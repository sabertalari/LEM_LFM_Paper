# Load Julia Packages
#--------------------
using PowerModels
using Ipopt
using JuMP

using NamedArrays
using LinearAlgebra
using Plots
plotlyjs()
using DataFrames

using XLSX
using JLD2
using FileIO

# load ADMM scripts
include("admm_scripts/ac_opf_single_branch_update.jl")
include("admm_scripts/ac_opf_single_bus_update.jl")
include("admm_scripts/ac_opf_lambda_update.jl")
include("admm_scripts/ac_opf_residual_update.jl")
include("admm_scripts/ac_opf_rho_update.jl")
include("admm_scripts/ac_opf_epsilon_update.jl")
include("ac_pf_central.jl")
include("ac_lfm_central.jl")
include("ac_lfm_admm.jl")

# Load System Data
# ----------------
powermodels_path = joinpath(dirname(pathof(PowerModels)), "..")

# Choose test case
test_case = "simple_mv_open_ring_net_lfm" # "simple_mv_closed_ring_net"  # 
# file_name = "$(powermodels_path)/test/data/matpower/$(test_case).m"
file_name = "./test_cases/$(test_case).m"
# note: change this string to modify the network data that will be loaded

# Load the data file
data = PowerModels.parse_file(file_name)

# Add zeros to turn linear objective functions into quadratic ones
# so that additional parameter checks are not required
PowerModels.standardize_cost_terms!(data, order=2)

# Adds reasonable rate_a values to branches without them
PowerModels.calc_thermal_limits!(data)

# Use build_ref to filter out inactive components
ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
# Note: ref contains all the relevant system parameters needed to build the OPF model
# When we introduce constraints and variable bounds below, we use the parameters in ref.

# Information from the LEM model
baseMVA_LEM = 2.3
number_ecs = 5
timestep = 16

ecm = Dict()
P_lem = Dict()
Q_lem = Dict()
R_pᵘ = Dict()
R_pᵈ = Dict()
R_qᵘ = Dict()
R_qᵈ = Dict()

# The first index is reserved for the generator at the slack bus
for ec in 2:(number_ecs+1)

    ecm[ec] = XLSX.readxlsx(abspath(joinpath(pwd(), "../Energy community management/ECM$((ec-1))-res.xlsx")))

    P_lem[ec] = ecm[ec]["P_EC!C2:C25"][timestep] * baseMVA_LEM / ref[:baseMVA]
    Q_lem[ec] = ecm[ec]["Q_EC!C2:C25"][timestep] * baseMVA_LEM / ref[:baseMVA]
    R_pᵘ[ec] = ecm[ec]["p_res_up!C2:C25"][timestep] * baseMVA_LEM / ref[:baseMVA]
    R_pᵈ[ec] = ecm[ec]["p_res_dn!C2:C25"][timestep] * baseMVA_LEM / ref[:baseMVA]
    R_qᵘ[ec] = ecm[ec]["q_res_up!C2:C25"][timestep] * baseMVA_LEM / ref[:baseMVA]
    R_qᵈ[ec] = ecm[ec]["q_res_dn!C2:C25"][timestep] * baseMVA_LEM / ref[:baseMVA]
end

# Read initial flex prices from LEM
ec_prices = XLSX.readxlsx(abspath(joinpath(pwd(), "../Energy community management/Prices in EC.xlsx")))
# print(ec_prices["Sheet1!E2:E25"]) # ECM3 + 4
# print(ec_prices["Sheet1!F2:F25"]) # ECM5
# print(ec_prices["Sheet1!G2:G25"]) # ECM1 + 2

ref[:gen][2]["cost"][2] = ec_prices["Sheet1!G2:G25"][timestep] * 1000
ref[:gen][3]["cost"][2] = ec_prices["Sheet1!G2:G25"][timestep] * 1000
ref[:gen][4]["cost"][2] = ec_prices["Sheet1!E2:E25"][timestep] * 1000
ref[:gen][5]["cost"][2] = ec_prices["Sheet1!E2:E25"][timestep] * 1000
ref[:gen][6]["cost"][2] = ec_prices["Sheet1!F2:F25"][timestep] * 1000

ref[:gen][1]["cost"][3] = 0.0
ref[:gen][2]["cost"][3] = 0.0
ref[:gen][3]["cost"][3] = 0.0
ref[:gen][4]["cost"][3] = 0.0
ref[:gen][5]["cost"][3] = 0.0
ref[:gen][6]["cost"][3] = 0.0

# Run powerflow with results from LEM to obtain power on the lines
# PowerModels.compute_ac_pf(data) # p_ij, q_ij = 
ac_pf_model, p_pf, q_pf, vm_pf, va_pf, pg_pf, qg_pf = central_ac_pf(ref, P_lem, Q_lem)
solution_summary(ac_pf_model, verbose=true)

# Add result of slackbus from powerflow to LEM data
P_lem[1] = pg_pf[1]
Q_lem[1] = qg_pf[1]


# Set flexibility demand from upstream network
fp_demand = -0.1 / ref[:baseMVA]
fq_demand = 0.0 / ref[:baseMVA]

# Solve centralized LFM for reference
fp_central, fq_central, fpl_central, fql_central, fvm_central, fva_central, cost_central, dual_p_central, dual_q_central, lfm_model_central= central_lfm(ref, P_lem, Q_lem, R_pᵘ, R_pᵈ, R_qᵘ, R_qᵈ, fp_demand, fq_demand, p_pf, q_pf, vm_pf, va_pf)
solution_summary(lfm_model_central, verbose=true)
if "$(termination_status(lfm_model_central))" == "LOCALLY_INFEASIBLE"
    throw(ErrorException("Central solution is LOCALLY_INFEASIBLE"))
end

K= 10000 # Max iterations

# Feasibility tolerances
ϵ_abs = 0.000001 # <= 10e-6 absolute tolerance for algorithm termination (Mhanna) | 0.001 (Gebbran)
ϵ_rel = 0.00002 # <= 5 * 10e-5 relative tolerance for algorithm termination (Mhanna) | 0.01 (Gebbran)

# Parameters for rho update
μ_decr = 1.15 # 1.15 (Gebbran)
μ_incr = 1/μ_decr # 1/μ_decr (Gebbran)
# μ_decr_lst = [1.13, 1.15, 1.17, 1.2]

τ_incr = 1.15 # 1.15 (Gebbran)
τ_decr = 0.9 # 0.9 (Gebbran)

# Upper and lower limit for penalty factor
ρ̅ = 10^5
ρ̲ = 5
# Iteration delay before updating ρ
if test_case == "simple_mv_closed_ring_net"
    ρ_delay = 10
elseif test_case == "simple_mv_open_ring_net_lfm"
    ρ_delay = 20
else
    throw(ErrorException("No ρ_delay defined for $(test_case)"))
end

@time begin
### Start of ADMM loop ###
fp, fq, f_plᵢ, f_qlᵢ, vm, va, f_pl, f_ql, vmₗ, vaₗ, r_norm, s_norm, ϵ_pri, ϵ_dual, cost, dual_p_admm, dual_q_admm, final_k, ρ = admm(ref, ϵ_abs, ϵ_rel, 
                                                                                                                                    μ_decr, μ_incr, 
                                                                                                                                    τ_incr, τ_decr, 
                                                                                                                                    ρ̅, ρ̲, ρ_delay, 
                                                                                                                                    fp_demand, fq_demand,
                                                                                                                                    p_pf, q_pf, vm_pf, va_pf,
                                                                                                                                    K)

                                                                                                                                
end

### ADMM performance ###
# Plot development of residuals and epsilons over iterations

df_flex_p = DataFrame(gen = Any[], Flexibility = Any[], dual = Any[])
df_flex_p_cent = DataFrame(gen = Any[], Flexibility = Any[], dual = Any[])
for (g, gen) in (ref[:gen])
    push!(df_flex_p,[ g, (fp[g, final_k]*ref[:baseMVA]), (dual_p_admm[ref[:gen][g]["gen_bus"], final_k]*ref[:baseMVA])])
    push!(df_flex_p_cent,[ g, (fp_central[g]*ref[:baseMVA]), (dual_p_central[g]*ref[:baseMVA])])
end

print("\nResults for ", ref[:name], ":\n")
print("ADMM:\n")
show(df_flex_p, allrows=true)
print("\nCentralized:\n")
show(df_flex_p_cent, allrows=true)
print("\n")

f_admm_sum = 0.0
f_cent_sum = 0.0
for g in keys(ref[:gen])
    if ref[:gen][g]["gen_bus"] in keys(ref[:ref_buses])
        continue
    else
        global f_admm_sum += fp[g, final_k]
        global f_cent_sum += fp_central[g]
    end
end
print("Flexibility Summations:\n")
print("ADMM: $(f_admm_sum) \nCentralized: $(f_cent_sum)")

plt_primal = plot(hcat(r_norm[1:final_k], ϵ_pri[1:final_k]),
                yaxis=(:log10, [0.00000001, :auto]),
                label=permutedims(["r", "ϵ_pri"]),
                legend=:outertopright,
                legendgroup="group1",
                legendgrouptitle_text="Primal Residual",
                fontfamily = "times",
                palette = :darktest)
plt_dual = plot(hcat(s_norm[1:final_k], ϵ_dual[1:final_k]),
                yaxis=(:log10, [0.00000001, :auto]),
                label=permutedims(["s", "ϵ_dual"]),
                legend=:outertopright,
                legendgroup="group2",
                legendgrouptitle_text="Dual Residual",
                fontfamily = "times",
                palette = :lightrainbow)
plt_flex = plot(transpose(dual_p_admm[3:7,1:final_k]/-1000), xlabel="Iterations", ylabel="€/kWh", 
                # yaxis=(:log10, [0.00000001, :auto]), 
                legend=:topright, 
                label=permutedims(["EC1","EC2","EC3","EC4","EC5"]),
                fontfamily = "times",
                legendgroup="group3",
                legendgrouptitle_text="Flexibility Prices")

plot(plt_primal,plt_dual,plt_flex, layout = grid(3, 1, heights=[0.2 ,0.2, 0.2]))
