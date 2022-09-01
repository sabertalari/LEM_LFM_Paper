function update_single_bus(ref, vmₗ, vaₗ, f_pl, f_ql, vmk, vak, fpk, fqk, f_plᵢk, f_qlᵢk, P_lem, Q_lem, R_pᵈ, R_pᵘ, R_qᵈ, R_qᵘ, λ_pl, λ_ql, λ_vm, λ_va, ρ_pl, ρ_ql, ρ_vm, ρ_va, i, p_pf, q_pf, fp_demand, fq_demand)

    # Initialize a JuMP Optimization Model
    #-------------------------------------
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    
    # model = Model(with_optimizer(Ipopt.Optimizer, linear_solver="ma27", print_level=0))
    # note: print_level changes the amount of solver information printed to the terminal

    # Add Optimization and State Variables
    # ------------------------------------

    # Add voltage angles va for each bus
    @variable(model, va)
    # note: [i in keys(ref[:bus])] adds one `va` variable for each bus in the network

    # Add voltage angles vm for each bus
    @variable(model, ref[:bus][i]["vmin"] <= vm <= ref[:bus][i]["vmax"]) # , start=1.0
    # @variable(model, vm[i in keys(ref[:bus])]) # , start=1.0
    # note: this vairable also includes the voltage magnitude limits and a starting value

    # Add flexibilities for the bus, limited by the results of the LEM model
    @variable(model, fp[g in ref[:bus_gens][i]])
    @variable(model, fq[g in ref[:bus_gens][i]])

    for g in ref[:bus_gens][i]
        if i in keys(ref[:ref_buses])
            @assert ref[:bus][i]["bus_type"] == 3
            @constraint(model, va == 0)
            @constraint(model, vm == 1)
            @constraint(model, fp[g] == fp_demand)
            @constraint(model, fq[g] == fq_demand)
        else
            @constraint(model, -R_pᵈ[g] <= fp[g] <= R_pᵘ[g])
            @constraint(model, -R_qᵈ[g] <= fq[g] <= R_qᵘ[g])
        end
    end
    

    # No more generation variables in LFM model, instead use LEM results
    pg = Dict()
    qg = Dict()
    for g in ref[:bus_gens][i]
        pg[g] = P_lem[g]
        qg[g] = Q_lem[g]
    end

    # Consensus Variables
    # Add power flow variables p to represent the active power flow for each branch
    @variable(model, f_plᵢ[(l,f,t) in ref[:bus_arcs][i]])
    # Add power flow variables q to represent the reactive power flow for each branch
    @variable(model, f_qlᵢ[(l,f,t) in ref[:bus_arcs][i]])

    # Set voltage starting values from previous iteration (first one from powerflow)
    set_start_value(vm, vmk)
    set_start_value(va, vak)

    for g in ref[:bus_gens][i]
        set_start_value(fp[g], fpk[g]) # set_start_value(pg[g], p_gk[g])
        set_start_value(fq[g], fqk[g]) # set_start_value(qg[g], q_gk[g])
    end
    
    for (l,f,t) in ref[:bus_arcs][i]
        set_start_value(f_plᵢ[(l,f,t)], f_plᵢk[string((l,f,t))])
        set_start_value(f_qlᵢ[(l,f,t)], f_qlᵢk[string((l,f,t))])
    end

    # Line Constraints
    for (l,f,t) in ref[:bus_arcs][i]
        @constraint(model, -ref[:branch][l]["rate_a"] <= (f_plᵢ[(l,f,t)] + p_pf[(l,f,t)]) <= ref[:branch][l]["rate_a"])
        @constraint(model, -ref[:branch][l]["rate_a"] <= (f_qlᵢ[(l,f,t)] + q_pf[(l,f,t)]) <= ref[:branch][l]["rate_a"])
    end

    # Nodal power balance constraints
    # Build a list of the loads and shunt elements connected to the bus i
    # bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
    bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

    @constraint(model, p_nodal_bal,
        sum(f_plᵢ[a] + p_pf[a] for a in ref[:bus_arcs][i]) ==   # sum of active power flow on lines from bus i + line usage from power flow calculation ==
        sum(pg[g] + fp[g] for g in ref[:bus_gens][i]) -         # sum of active power generation at bus i + possible active flexibility at the bus -
        # sum(load["pd"] for load in bus_loads) -                 # sum of active load consumption at bus i -
        sum(shunt["gs"] for shunt in bus_shunts)*vm^2        # sum of active shunt element injections at bus i
    )

    # Reactive power balance at node i
    @constraint(model, q_nodal_bal,
        sum(f_qlᵢ[a] + q_pf[a] for a in ref[:bus_arcs][i]) ==     # sum of reactive power flow on lines from bus i + line usage from power flow calculation ==
        sum(qg[g] + fq[g] for g in ref[:bus_gens][i]) +                 # sum of reactive power generation at bus i - possible reactive flexibility at bus i -
        # sum(load["qd"] for load in bus_loads) +                 # sum of reactive load consumption at bus i +
        sum(shunt["bs"] for shunt in bus_shunts)*vm^2        # sum of reactive shunt element injections at bus i
    )

    # Add Objective Function
    # ----------------------

    # Minimize the cost of active power generation
    # assumes costs are given as quadratic functions
    # pl, ql, val, vml from branch update

    # plᵢ - pl
    @objective(model, Min,
        (sum(ref[:gen][g]["cost"][1]*fp[g]^2 + ref[:gen][g]["cost"][2]*fp[g] + ref[:gen][g]["cost"][3] for g in ref[:bus_gens][i])
        + sum(λ_pl[string(a)] * (f_plᵢ[a] - f_pl[string(a)])
            + λ_ql[string(a)] * (f_qlᵢ[a] - f_ql[string(a)])
            + λ_vm[string(a)] * (vm - vmₗ[string(a)])
            + λ_va[string(a)] * (va - vaₗ[string(a)]) 
            for a in ref[:bus_arcs][i])
        + (
            sum(ρ_pl[string(a)]/2 * (f_plᵢ[a] - f_pl[string(a)])^2
                + ρ_ql[string(a)]/2 * (f_qlᵢ[a] - f_ql[string(a)])^2
                + ρ_vm[string(a)]/2 * (vm - vmₗ[string(a)])^2
                + ρ_va[string(a)]/2 * (va - vaₗ[string(a)])^2
                for a in ref[:bus_arcs][i])
            )
        )
    )

    optimize!(model)


    if length(ref[:bus_gens][i]) > 0
        cost = sum(ref[:gen][g]["cost"][1]*JuMP.value.(fp)[g]^2 + ref[:gen][g]["cost"][2]*JuMP.value.(fp)[g] + ref[:gen][g]["cost"][3] for g in ref[:bus_gens][i])
    else
        cost = 0.0
    end

    dual_p = dual(p_nodal_bal)
    dual_q = dual(q_nodal_bal)

    return JuMP.value.(fp), JuMP.value.(fq), JuMP.value.(f_plᵢ), JuMP.value.(f_qlᵢ), JuMP.value.(vm), JuMP.value.(va), dual_p, dual_q, cost
end