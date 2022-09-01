function update_single_branch(ref, vm, va, f_plᵢ, f_qlᵢ, vmₗk, vaₗk, f_plk, f_qlk, p_pf, q_pf, λ_pl, λ_ql, λ_vm, λ_va, ρ_pl, ρ_ql, ρ_vm, ρ_va, l, i, j)

    # Initialize a JuMP Optimization Model
    #-------------------------------------
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    
    # model = Model(with_optimizer(Ipopt.Optimizer, linear_solver="ma27", print_level=0))
    # note: print_level changes the amount of solver information printed to the terminal

    # Add power flow variables f_pl to represent the active power flow for each branch due to flexibility allocation
    @variable(model, f_pl[a in [(l,i,j),(l,j,i)]])
    # Add power flow variables f_ql to represent the reactive power flow for each branch due to flexibility allocation
    @variable(model, f_ql[a in [(l,i,j),(l,j,i)]])
    # note: ref[:arcs] includes both the from (i,j) and the to (j,i) sides of a branch

    # Consensus Variables
    # Add voltage angles va for each bus
    @variable(model, vaₗ[a in [(l,i,j),(l,j,i)]])
    # note: [i in keys(ref[:bus])] adds one `va` variable for each bus in the network

    # Add voltage angles vm for each bus
    @variable(model, ref[:bus][i]["vmin"] <= vmₗ[a in [(l,i,j),(l,j,i)]] <= ref[:bus][i]["vmax"])  # , start=1.0
    # note: this vairable also includes the voltage magnitude limits and a starting value

    for a in [(l,i,j),(l,j,i)]
        set_start_value(vmₗ[a], vmₗk[string(a)])
        set_start_value(vaₗ[a], vaₗk[string(a)])

        set_start_value(f_pl[a], f_plk[string(a)])
        set_start_value(f_ql[a], f_qlk[string(a)])
    end

    # Fix the voltage angle to zero (and magnitude to one) at the branch,
    # connectd to the reference bus
    for (b,bus) in ref[:ref_buses]
        # Check if either side of the current branch is connected to the reference bus
        if (l,i,j) in ref[:bus_arcs][b]
            @constraint(model, vaₗ[(l,i,j)] == 0)
            @constraint(model, vmₗ[(l,i,j)] == 1)  # ?
        elseif (l,j,i) in ref[:bus_arcs][b]
            @constraint(model, vaₗ[(l,j,i)] == 0)
            @constraint(model, vmₗ[(l,j,i)] == 1)  # ?
        end
    end

    # Powerflow constraints, usually defined with the variable p_pf. Not possible due to addition
    for (l,f,t) in [(l,i,j),(l,j,i)]
        @constraint(model, -ref[:branch][l]["rate_a"] <= (f_pl[(l,f,t)] + p_pf[(l,f,t)]) <= ref[:branch][l]["rate_a"])
        @constraint(model, -ref[:branch][l]["rate_a"] <= (f_ql[(l,f,t)] + q_pf[(l,f,t)]) <= ref[:branch][l]["rate_a"])
    end

    # Branch power flow physics and limit constraints
    # Build the from variable id of the i-th branch, which is a tuple given by (branch id, from bus, to bus)
    f_idx = (l,i,j) # (i, branch["f_bus"], branch["t_bus"])
    # Build the to variable id of the i-th branch, which is a tuple given by (branch id, to bus, from bus)
    t_idx = (l,j,i) # (i, branch["t_bus"], branch["f_bus"])
    # note: it is necessary to distinguish between the from and to sides of a branch due to power losses

    # Compute the branch parameters and transformer ratios from the data
    g, b = PowerModels.calc_branch_y(ref[:branch][l]) # PowerModels.calc_branch_y(branch)
    tr, ti = PowerModels.calc_branch_t(ref[:branch][l]) # PowerModels.calc_branch_t(branch)
    g_fr = ref[:branch][l]["g_fr"] # branch["g_fr"]
    b_fr = ref[:branch][l]["b_fr"] # branch["b_fr"]
    g_to = ref[:branch][l]["g_to"] # branch["g_to"]
    b_to = ref[:branch][l]["b_to"] # branch["b_to"]
    tm = ref[:branch][l]["tap"]^2 # branch["tap"]^2
    # note: tap is assumed to be 1.0 on non-transformer branches

    # AC Power Flow Constraints
    # From side of the branch flow (summation)
    @NLconstraint(model, f_pl[f_idx] + p_pf[f_idx] ==  (g+g_fr)/tm*vmₗ[f_idx]^2 + (-g*tr+b*ti)/tm*(vmₗ[f_idx]*vmₗ[t_idx]*cos(vaₗ[f_idx]-vaₗ[t_idx])) + (-b*tr-g*ti)/tm*(vmₗ[f_idx]*vmₗ[t_idx]*sin(vaₗ[f_idx]-vaₗ[t_idx])))
    @NLconstraint(model, f_ql[f_idx] + q_pf[f_idx]  == -(b+b_fr)/tm*vmₗ[f_idx]^2 - (-b*tr-g*ti)/tm*(vmₗ[f_idx]*vmₗ[t_idx]*cos(vaₗ[f_idx]-vaₗ[t_idx])) + (-g*tr+b*ti)/tm*(vmₗ[f_idx]*vmₗ[t_idx]*sin(vaₗ[f_idx]-vaₗ[t_idx])))

    # To side of the branch flow (summation)
    @NLconstraint(model, f_pl[t_idx] + p_pf[t_idx] == (g+g_to)*vmₗ[t_idx]^2 + (-g*tr-b*ti)/tm*(vmₗ[t_idx]*vmₗ[f_idx]*cos(vaₗ[t_idx]-vaₗ[f_idx])) + (-b*tr+g*ti)/tm*(vmₗ[t_idx]*vmₗ[f_idx]*sin(vaₗ[t_idx]-vaₗ[f_idx])) )
    @NLconstraint(model, f_ql[t_idx] + q_pf[t_idx] == -(b+b_to)*vmₗ[t_idx]^2 - (-b*tr+g*ti)/tm*(vmₗ[t_idx]*vmₗ[f_idx]*cos(vaₗ[f_idx]-vaₗ[t_idx])) + (-g*tr-b*ti)/tm*(vmₗ[t_idx]*vmₗ[f_idx]*sin(vaₗ[t_idx]-vaₗ[f_idx])) )

    # Voltage angle difference limit
    @constraint(model, vaₗ[f_idx] - vaₗ[t_idx] <= ref[:branch][l]["angmax"])
    @constraint(model, vaₗ[f_idx] - vaₗ[t_idx] >= ref[:branch][l]["angmin"])

    # Apparent power limit, from side and to side
    @constraint(model, (f_pl[f_idx] + p_pf[f_idx])^2 + (f_ql[f_idx] + q_pf[f_idx])^2 <= ref[:branch][l]["rate_a"]^2)
    @constraint(model, (f_pl[t_idx] + p_pf[t_idx])^2 + (f_ql[t_idx] + q_pf[t_idx])^2 <= ref[:branch][l]["rate_a"]^2)

    @objective(model, Min,
        (sum(λ_pl[string((b,f,t))] * (f_plᵢ[string((b,f,t))] - f_pl[(b,f,t)])
            + λ_ql[string((b,f,t))] * (f_qlᵢ[string((b,f,t))] - f_ql[(b,f,t)])
            + λ_vm[string((b,f,t))] * (vm[string(f)] - vmₗ[(b,f,t)])
            + λ_va[string((b,f,t))] * (va[string(f)] - vaₗ[(b,f,t)]) 
            for (b,f,t) in [(l,i,j),(l,j,i)])
        + (
            sum(ρ_pl[string((b,f,t))]/2 * (f_plᵢ[string((b,f,t))] - f_pl[(b,f,t)])^2
                + ρ_ql[string((b,f,t))]/2 * (f_qlᵢ[string((b,f,t))] - f_ql[(b,f,t)])^2
                + ρ_vm[string((b,f,t))]/2 * (vm[string(f)] - vmₗ[(b,f,t)])^2
                + ρ_va[string((b,f,t))]/2 * (va[string(f)] - vaₗ[(b,f,t)])^2
                for (b,f,t) in [(l,i,j),(l,j,i)])
            )
        )
    )
    # Solve the optimization problem
    optimize!(model)

    return JuMP.value.(f_pl), JuMP.value.(f_ql), JuMP.value.(vmₗ), JuMP.value.(vaₗ)
end
