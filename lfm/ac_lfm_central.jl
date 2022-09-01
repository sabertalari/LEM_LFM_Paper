function central_lfm(ref, P_lem, Q_lem, R_pᵘ, R_pᵈ, R_qᵘ, R_qᵈ, fp_demand, fq_demand, p_pf, q_pf, vm_pf, va_pf)
    # Initialize a JuMP Optimization Model
    #-------------------------------------
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)

    # Add voltage angles va for each bus
    @variable(model, va[i in keys(ref[:bus])])
    # note: [i in keys(ref[:bus])] adds one `va` variable for each bus in the network

    # Add voltage magnitude vm for each bus
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"])
    # note: this vairable also includes the voltage magnitude limits

    # Add flexibilities for the bus
    @variable(model, fp[g in keys(ref[:gen])])
    @variable(model, fq[g in keys(ref[:gen])])

    # Add power flow variables f_pl to represent the active power flow for each branch due to flexibility allocation
    @variable(model, f_pl[(l,i,j) in ref[:arcs]])
    # Add power flow variables f_ql to represent the reactive power flow for each branch due to flexibility allocation
    @variable(model, f_ql[(l,i,j) in ref[:arcs]])

    # No more generation variables in LFM model, instead use LEM results
    pg = Dict()
    qg = Dict()
    for g in keys(ref[:gen])
        pg[g] = P_lem[g]
        qg[g] = Q_lem[g]
    end

    # Limit flexibility of ECM buses by the results of the LEM model
    for g in keys(ref[:gen])
        if ref[:gen][g]["gen_bus"] in keys(ref[:ref_buses])
            @assert ref[:bus][ref[:gen][g]["gen_bus"]]["bus_type"] == 3
            continue
        else
            @constraint(model, -R_pᵈ[g] <= fp[g] <= R_pᵘ[g])
            @constraint(model, -R_qᵈ[g] <= fq[g] <= R_qᵘ[g])
        end
    end
    for i in keys(ref[:bus])
        set_start_value(vm[i], vm_pf[i])
        set_start_value(va[i], va_pf[i])
    end
    # Fix the voltage angle and magnitude at the reference bus
    for i in keys(ref[:ref_buses])
        @assert ref[:bus][i]["bus_type"] == 3
        @constraint(model, va[i] == 0)
        @constraint(model, vm[i] == 1)  # ?
        for g in ref[:bus_gens][i]
            @constraint(model, fp[g] == fp_demand)
            @constraint(model, fq[g] == fq_demand)
        end
    end


    p_nodal_bal = Dict()
    q_nodal_bal = Dict()
    # Nodal power balance constraints
    # Build a list of the loads and shunt elements connected to the bus i
    for i in keys(ref[:bus])
        # bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]


        # Active power balance at node i
        p_nodal_bal[i] = @constraint(model,
            sum(f_pl[a] + p_pf[a] for a in ref[:bus_arcs][i]) ==   # sum of active power flow on lines from bus i + line usage from power flow calculation ==
            sum(pg[g] + fp[g] for g in ref[:bus_gens][i]) -         # sum of active power generation at bus i + possible active flexibility at the bus -
            # sum(load["pd"] for load in bus_loads) -                 # sum of active load consumption at bus i -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2        # sum of active shunt element injections at bus i
        )

        # Reactive power balance at node i
        q_nodal_bal[i] = @constraint(model,
            sum(f_ql[a] + q_pf[a] for a in ref[:bus_arcs][i]) ==     # sum of reactive power flow on lines from bus i + line usage from power flow calculation ==
            sum(qg[g] + fq[g] for g in ref[:bus_gens][i]) +                 # sum of reactive power generation at bus i - possible reactive flexibility at bus i -  ## + because no load considered!
            # sum(load["qd"] for load in bus_loads) +                 # sum of reactive load consumption at bus i +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2        # sum of reactive shunt element injections at bus i
        )
    end

    # Powerflow constraints, usually defined with the variable p_pf. Not possible due to addition
    for (l,i,j) in ref[:arcs]
        @constraint(model, -ref[:branch][l]["rate_a"] <= (f_pl[(l,i,j)] + p_pf[(l,i,j)]) <= ref[:branch][l]["rate_a"])
        @constraint(model, -ref[:branch][l]["rate_a"] <= (f_ql[(l,i,j)] + q_pf[(l,i,j)]) <= ref[:branch][l]["rate_a"])
    end

    # Branch power flow physics and limit constraints
    for (i, branch) in ref[:branch]
        # Build the from variable id of the i-th branch, which is a tuple given by (branch id, from bus, to bus)
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        # Build the to variable id of the i-th branch, which is a tuple given by (branch id, to bus, from bus)
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        # note: it is necessary to distinguish between the from and to sides of a branch due to power losses

        # Compute the branch parameters and transformer ratios from the data
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2
        # note: tap is assumed to be 1.0 on non-transformer branches

        vm_fr = vm[branch["f_bus"]]         # vm_fr is a reference to the optimization variable vm on the from side of the branch
        vm_to = vm[branch["t_bus"]]         # vm_to is a reference to the optimization variable vm on the to side of the branch
        va_fr = va[branch["f_bus"]]         # va_fr is a reference to the optimization variable va on the from side of the branch
        va_to = va[branch["t_bus"]]         # va_to is a reference to the optimization variable va on the to side of the branch

        # AC Power Flow Constraints
        # From side of the branch flow (summation)
        @NLconstraint(model, f_pl[f_idx] + p_pf[f_idx] ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to))) 
        @NLconstraint(model, f_ql[f_idx] + q_pf[f_idx]  == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)))

        # To side of the branch flow (summation)
        @NLconstraint(model, f_pl[t_idx] + p_pf[t_idx] ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)))
        @NLconstraint(model, f_ql[t_idx] + q_pf[t_idx] == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)))

        # Voltage angle difference limit
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent power limit, from side and to side
        @constraint(model, (f_pl[f_idx] + p_pf[f_idx])^2 + (f_ql[f_idx] + q_pf[f_idx])^2 <= branch["rate_a"]^2)
        @constraint(model, (f_pl[t_idx] + p_pf[t_idx])^2 + (f_ql[t_idx] + q_pf[t_idx])^2 <= branch["rate_a"]^2)
    end

    # Add Objective Function
    # ----------------------

    # Minimize the cost of active power generation
    # assumes costs are given as quadratic functions

    @objective(model, Min,
        sum(ref[:gen][g]["cost"][1]*fp[g]^2 + ref[:gen][g]["cost"][2]*fp[g] + ref[:gen][g]["cost"][3] for g in keys(ref[:gen]))

    )

    optimize!(model)
    @info("Central LFM terminated with status: $(termination_status(model))")

    cost = sum(ref[:gen][g]["cost"][1]*JuMP.value.(fp)[g]^2 + ref[:gen][g]["cost"][2]*JuMP.value.(fp)[g] + ref[:gen][g]["cost"][3] for g in keys(ref[:gen]))
    dual_p = Dict()
    dual_q = Dict()
    fp_central = Dict()
    fq_central = Dict()
    fpl_central = Dict()
    fql_central = Dict()
    fvm_central = Dict()
    fva_central = Dict()
    for i in keys(ref[:bus])
        dual_p[i] =  dual(p_nodal_bal[i])
        dual_q[i] =  dual(q_nodal_bal[i])
        fvm_central[i] = JuMP.value(vm[i])
        fva_central[i] = JuMP.value(va[i])
        for g in ref[:bus_gens][i]
            fp_central[g] = JuMP.value(fp[g])
            fq_central[g] = JuMP.value(fq[g])
        end

        for a in ref[:bus_arcs][i]
            fpl_central[a] = JuMP.value(f_pl[a])
            fql_central[a] = JuMP.value(f_ql[a])
        end
    end
   


    return fp_central, fq_central, fpl_central, fql_central, fvm_central, fva_central, cost, dual_p, dual_q, model
end