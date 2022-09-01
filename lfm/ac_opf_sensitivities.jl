# Load Julia Packages
#--------------------
using PowerModels
using Ipopt
using JuMP

using NamedArrays
using LinearAlgebra
using Plots
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
test_case = "simple_mv_closed_ring_net" # "simple_mv_open_ring_net_lfm" # 
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
fp_central, fq_central, fpl_central, fql_central, fvm_central, fva_central, fcost_central, dual_p_central, dual_q_central, lfm_model_central= central_lfm(ref, P_lem, Q_lem, R_pᵘ, R_pᵈ, R_qᵘ, R_qᵈ, fp_demand, fq_demand, p_pf, q_pf, vm_pf, va_pf)
solution_summary(lfm_model_central, verbose=true)

K= 10000 # Max iterations

# Feasibility tolerances
ϵ_abs = 0.000001 # <= 10e-6 absolute tolerance for algorithm termination (Mhanna) | 0.001 (Gebbran)
ϵ_rel = 0.00002 # <= 5 * 10e-5 relative tolerance for algorithm termination (Mhanna) | 0.01 (Gebbran)
ϵ_abs_lst = [0.001, 0.0001, 0.00001, 0.000001,  0.0000001]
ϵ_rel_lst = [0.02,  0.002,  0.0002,  0.00002,   0.000002]

# Parameters for rho update
μ_decr = 1.15 # 1.15 (Gebbran)
μ_incr = 1/μ_decr # 1/μ_decr (Gebbran)
# μ_decr_lst = [1.13, 1.15, 1.17, 1.2]

τ_incr = 1.15 # 1.15 (Gebbran)
τ_decr = 0.9 # 0.9 (Gebbran)

τ_incr_lst = [1.05, 1.1, 1.15, 1.2, 1.25, 1.3] # 1.15 (Gebbran)
τ_decr_lst = [0.8, 0.85, 0.9, 0.95, 0.1, 1.05] # 0.9 (Gebbran)

# Upper and lower limit for penalty factor
ρ̅ = 10^5
ρ̲ = 5
# Iteration delay before updating ρ
ρ_delay = 30
ρ_delay_lst = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 150]

### Start of ADMM loop ###

ρ_results = Dict()
for ρ_delay in ρ_delay_lst

    print("\n")
    @info("Started ADMM with ρ_delay: $(ρ_delay)")
    fp, fq, f_plᵢ, f_qlᵢ, vm, va, f_pl, f_ql, vmₗ, vaₗ, r_norm, s_norm, ϵ_pri, ϵ_dual, cost, dual_p_admm, dual_q_admm, final_k, ρ = admm(ref, ϵ_abs, ϵ_rel, 
                                                                                                                                        μ_decr, μ_incr, 
                                                                                                                                        τ_incr, τ_decr, 
                                                                                                                                        ρ̅, ρ̲, ρ_delay, 
                                                                                                                                        fp_demand, fq_demand,
                                                                                                                                        p_pf, q_pf, vm_pf, va_pf,
                                                                                                                                        K)

    ρ_results["ρ_delay-$(ρ_delay)"] = Dict()
    ρ_results["ρ_delay-$(ρ_delay)"]["fp"] = fp
    ρ_results["ρ_delay-$(ρ_delay)"]["fq"] = fq
    ρ_results["ρ_delay-$(ρ_delay)"]["f_plᵢ"] = f_plᵢ
    ρ_results["ρ_delay-$(ρ_delay)"]["f_qlᵢ"] = f_qlᵢ
    ρ_results["ρ_delay-$(ρ_delay)"]["vm"] = vm
    ρ_results["ρ_delay-$(ρ_delay)"]["va"] = va
    ρ_results["ρ_delay-$(ρ_delay)"]["f_pl"] = f_pl
    ρ_results["ρ_delay-$(ρ_delay)"]["f_ql"] = f_ql
    ρ_results["ρ_delay-$(ρ_delay)"]["vmₗ"] = vmₗ
    ρ_results["ρ_delay-$(ρ_delay)"]["vaₗ"] = vaₗ
    ρ_results["ρ_delay-$(ρ_delay)"]["r_norm"] = r_norm
    ρ_results["ρ_delay-$(ρ_delay)"]["s_norm"] = s_norm
    ρ_results["ρ_delay-$(ρ_delay)"]["ϵ_pri"] = ϵ_pri
    ρ_results["ρ_delay-$(ρ_delay)"]["ϵ_dual"] = ϵ_dual
    ρ_results["ρ_delay-$(ρ_delay)"]["cost"] = cost
    ρ_results["ρ_delay-$(ρ_delay)"]["final_k"] = final_k
    ρ_results["ρ_delay-$(ρ_delay)"]["ρ"] = ρ

    print("Final iter.: $(final_k) > ρ $(ρ[final_k]) r_norm $(r_norm[final_k]) ϵ_pri $(ϵ_pri[final_k])| s_norm $(s_norm[final_k]) ϵ_dual $(ϵ_dual[final_k]) | cost $(sum(cost[:,final_k]))")
end

ϵ_results = Dict()
for ϵ_nr in eachindex(ϵ_abs_lst)

    ϵ_abs = ϵ_abs_lst[ϵ_nr]
    ϵ_rel = ϵ_rel_lst[ϵ_nr]

    print("\n")
    @info("Started ADMM with ϵ_abs=$(ϵ_abs) and ϵ_rel=$(ϵ_rel)")
    fp, fq, f_plᵢ, f_qlᵢ, vm, va, f_pl, f_ql, vmₗ, vaₗ, r_norm, s_norm, ϵ_pri, ϵ_dual, cost, dual_p_admm, dual_q_admm, final_k, ρ = admm(ref, ϵ_abs, ϵ_rel, 
                                                                                                                                        μ_decr, μ_incr, 
                                                                                                                                        τ_incr, τ_decr, 
                                                                                                                                        ρ̅, ρ̲, ρ_delay, 
                                                                                                                                        fp_demand, fq_demand,
                                                                                                                                        p_pf, q_pf, vm_pf, va_pf,
                                                                                                                                        K)

    ϵ_results["ϵ_abs-$(ϵ_abs)"] = Dict()
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["ϵ_abs"] = ϵ_abs
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["ϵ_rel"] = ϵ_rel
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["fp"] = fp
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["fq"] = fq
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["f_plᵢ"] = f_plᵢ
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["f_qlᵢ"] = f_qlᵢ
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["vm"] = vm
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["va"] = va
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["f_pl"] = f_pl
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["f_ql"] = f_ql
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["vmₗ"] = vmₗ
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["vaₗ"] = vaₗ
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["r_norm"] = r_norm
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["s_norm"] = s_norm
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["ϵ_pri"] = ϵ_pri
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["ϵ_dual"] = ϵ_dual
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["cost"] = cost
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["final_k"] = final_k
    ϵ_results["ϵ_abs-$(ϵ_abs)"]["ρ"] = ρ

    print("Final iter.: $(final_k) > ρ $(ρ[final_k]) r_norm $(r_norm[final_k]) ϵ_pri $(ϵ_pri[final_k])| s_norm $(s_norm[final_k]) ϵ_dual $(ϵ_dual[final_k]) | cost $(sum(cost[:,final_k]))")
end

τ_results = Dict()
for τ_nr in eachindex(τ_incr_lst)

    τ_incr = τ_incr_lst[τ_nr]
    τ_decr = τ_decr_lst[τ_nr]

    print("\n")
    @info("Started ADMM with τ_incr=$(τ_incr) and τ_decr=$(τ_decr)")
    fp, fq, f_plᵢ, f_qlᵢ, vm, va, f_pl, f_ql, vmₗ, vaₗ, r_norm, s_norm, ϵ_pri, ϵ_dual, cost, dual_p_admm, dual_q_admm, final_k, ρ = admm(ref, ϵ_abs, ϵ_rel, 
                                                                                                                                        μ_decr, μ_incr, 
                                                                                                                                        τ_incr, τ_decr, 
                                                                                                                                        ρ̅, ρ̲, ρ_delay, 
                                                                                                                                        fp_demand, fq_demand,
                                                                                                                                        p_pf, q_pf, vm_pf, va_pf,
                                                                                                                                        K)

    τ_results["τ_incr-$(τ_incr)"] = Dict()
    τ_results["τ_incr-$(τ_incr)"]["τ_incr"] = τ_incr
    τ_results["τ_incr-$(τ_incr)"]["τ_decr"] = τ_decr
    τ_results["τ_incr-$(τ_incr)"]["fp"] = fp
    τ_results["τ_incr-$(τ_incr)"]["fq"] = fq
    τ_results["τ_incr-$(τ_incr)"]["f_plᵢ"] = f_plᵢ
    τ_results["τ_incr-$(τ_incr)"]["f_qlᵢ"] = f_qlᵢ
    τ_results["τ_incr-$(τ_incr)"]["vm"] = vm
    τ_results["τ_incr-$(τ_incr)"]["va"] = va
    τ_results["τ_incr-$(τ_incr)"]["f_pl"] = f_pl
    τ_results["τ_incr-$(τ_incr)"]["f_ql"] = f_ql
    τ_results["τ_incr-$(τ_incr)"]["vmₗ"] = vmₗ
    τ_results["τ_incr-$(τ_incr)"]["vaₗ"] = vaₗ
    τ_results["τ_incr-$(τ_incr)"]["r_norm"] = r_norm
    τ_results["τ_incr-$(τ_incr)"]["s_norm"] = s_norm
    τ_results["τ_incr-$(τ_incr)"]["ϵ_pri"] = ϵ_pri
    τ_results["τ_incr-$(τ_incr)"]["ϵ_dual"] = ϵ_dual
    τ_results["τ_incr-$(τ_incr)"]["cost"] = cost
    τ_results["τ_incr-$(τ_incr)"]["final_k"] = final_k
    τ_results["τ_incr-$(τ_incr)"]["ρ"] = ρ

    print("Final iter.: $(final_k) > ρ $(ρ[final_k]) r_norm $(r_norm[final_k]) ϵ_pri $(ϵ_pri[final_k])| s_norm $(s_norm[final_k]) ϵ_dual $(ϵ_dual[final_k]) | cost $(sum(cost[:,final_k]))")
end


### ADMM performance ###
# Plot development of residuals and epsilons over iterations
if length(ρ_delay_lst) > 0
    for ρ_delay in ρ_delay_lst
        final_k = ρ_results["ρ_delay-$(ρ_delay)"]["final_k"]
        df_res = DataFrame(r = ρ_results["ρ_delay-$(ρ_delay)"]["r_norm"][1:final_k], ϵ_pri = ρ_results["ρ_delay-$(ρ_delay)"]["ϵ_pri"][1:final_k], 
                            s = ρ_results["ρ_delay-$(ρ_delay)"]["s_norm"][1:final_k], ϵ_dual = ρ_results["ρ_delay-$(ρ_delay)"]["ϵ_dual"][1:final_k])
        plt = plot(Matrix(df_res), labels=permutedims(names(df_res)), legend=:outertopright, yaxis=(:log10, [0.00000001, :auto]), layout=(4,1))
        plt[:plot_title] = " ρ_delay: $(ρ_delay)" # ref[:name] * 
        display(plot(plt))
    end

    df_ρ = DataFrame(ρ_delay = Int[], Iterations = Int[], Gap= Any[])
    for ρ_delay in ρ_delay_lst

        final_k = ρ_results["ρ_delay-$(ρ_delay)"]["final_k"]
        f_admm_sum = 0.0
        f_cent_sum = 0.0
        for g in keys(ref[:gen])
            if ref[:gen][g]["gen_bus"] in keys(ref[:ref_buses])
                continue
            else
                f_admm_sum += ρ_results["ρ_delay-$(ρ_delay)"]["fp"][g, final_k+1]
                f_cent_sum += fp_central[g]
            end
        end
        push!(df_ρ,[ρ_delay, final_k, ((f_cent_sum-f_admm_sum)/f_cent_sum*100)])
    end
    show(df_ρ, allrows=true)
end

if length(ϵ_abs_lst) > 0
    for ϵ_abs in ϵ_abs_lst
        final_k = ϵ_results["ϵ_abs-$(ϵ_abs)"]["final_k"]
        ϵ_rel = ϵ_results["ϵ_abs-$(ϵ_abs)"]["ϵ_rel"]
        df_res = DataFrame(r = ϵ_results["ϵ_abs-$(ϵ_abs)"]["r_norm"][1:final_k], ϵ_pri = ϵ_results["ϵ_abs-$(ϵ_abs)"]["ϵ_pri"][1:final_k], 
                            s = ϵ_results["ϵ_abs-$(ϵ_abs)"]["s_norm"][1:final_k], ϵ_dual = ϵ_results["ϵ_abs-$(ϵ_abs)"]["ϵ_dual"][1:final_k])
        plt = plot(Matrix(df_res), labels=permutedims(names(df_res)), legend=:outertopright, yaxis=(:log10, [0.00000001, :auto]), layout=(4,1))
        plt[:plot_title] = " ϵ_abs-$(ϵ_abs) / ϵ_rel=$(ϵ_rel)" # ref[:name] * 
        display(plot(plt))
    end

    df_ϵ = DataFrame(ϵ_abs = Any[], ϵ_rel = Any[], Iterations = Int[], Gap= Any[])
    for ϵ_abs in ϵ_abs_lst

        ϵ_rel = ϵ_results["ϵ_abs-$(ϵ_abs)"]["ϵ_rel"]
        final_k = ϵ_results["ϵ_abs-$(ϵ_abs)"]["final_k"]
        f_admm_sum = 0.0
        f_cent_sum = 0.0
        for g in keys(ref[:gen])
            if ref[:gen][g]["gen_bus"] in keys(ref[:ref_buses])
                continue
            else
                f_admm_sum += ϵ_results["ϵ_abs-$(ϵ_abs)"]["fp"][g, final_k+1]
                f_cent_sum += fp_central[g]
            end
        end
        push!(df_ϵ,[ϵ_abs, ϵ_rel, final_k, ((f_cent_sum-f_admm_sum)/f_cent_sum*100)])
    end
    show(df_ϵ, allrows=true)
end

if length(τ_incr_lst) > 0
    for τ_incr in τ_incr_lst
        final_k = τ_results["τ_incr-$(τ_incr)"]["final_k"]
        τ_decr = τ_results["τ_incr-$(τ_incr)"]["τ_decr"]
        df_res = DataFrame(r = τ_results["τ_incr-$(τ_incr)"]["r_norm"][1:final_k], ϵ_pri = τ_results["τ_incr-$(τ_incr)"]["ϵ_pri"][1:final_k], 
                            s = τ_results["τ_incr-$(τ_incr)"]["s_norm"][1:final_k], ϵ_dual = τ_results["τ_incr-$(τ_incr)"]["ϵ_dual"][1:final_k])
        plt = plot(Matrix(df_res), labels=permutedims(names(df_res)), legend=:outertopright, yaxis=(:log10, [0.00000001, :auto]), layout=(4,1))
        plt[:plot_title] = " τ_incr-$(τ_incr) / τ_decr=$(τ_decr)" # ref[:name] * 
        display(plot(plt))
    end

    df_τ = DataFrame(τ_incr = Any[], τ_decr = Any[], Iterations = Int[], Gap= Any[])
    for τ_incr in τ_incr_lst

        τ_decr = τ_results["τ_incr-$(τ_incr)"]["τ_decr"]
        final_k = τ_results["τ_incr-$(τ_incr)"]["final_k"]
        f_admm_sum = 0.0
        f_cent_sum = 0.0
        for g in keys(ref[:gen])
            if ref[:gen][g]["gen_bus"] in keys(ref[:ref_buses])
                continue
            else
                f_admm_sum += τ_results["τ_incr-$(τ_incr)"]["fp"][g, final_k+1]
                f_cent_sum += fp_central[g]
            end
        end
        push!(df_τ,[τ_incr, τ_decr, final_k, ((f_cent_sum-f_admm_sum)/f_cent_sum*100)])
    end
    show(df_τ, allrows=true)
end
