"""
Given a JuMP model and a PowerModels network data structure,
Builds an AC-PF formulation of the given data and returns the JuMP model
"""
function central_ac_pf(ref, P_lem, Q_lem, model=Model(with_optimizer(Ipopt.Optimizer, print_level=0)))

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, vm[i in keys(ref[:bus])], start=1.0)

    @variable(model, pg[i in keys(ref[:gen])])
    @variable(model, qg[i in keys(ref[:gen])])

    for i in keys(ref[:gen])
        if ref[:gen][i]["gen_bus"] in keys(ref[:ref_buses])
            @assert ref[:bus][i]["bus_type"] == 3
            continue
        else
            @constraint(model, pg[i] == P_lem[i])
            @constraint(model, qg[i] == Q_lem[i])
        end
    end

    p = Dict()
    q = Dict()
    for (i,branch) in ref[:branch]
        # AC Line Flow expressions
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        p[f_idx] = @NLexpression(model,  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        q[f_idx] = @NLexpression(model, -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        p[t_idx] = @NLexpression(model,  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        q[t_idx] = @NLexpression(model, -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
    end


    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @assert bus["bus_type"] == 3
        @constraint(model, va[i] == 0)
        @constraint(model, vm[i] == bus["vm"])
    end

    for (i,bus) in ref[:bus]
        # bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        @NLconstraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            # sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )

        @NLconstraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) +
            # sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )

        # PV Bus Constraints
        # if length(ref[:bus_gens][i]) > 0 && !(i in keys(ref[:ref_buses]))
        #     # this assumes inactive generators are filtered out of bus_gens
        #     @assert bus["bus_type"] == 2

        #     @constraint(model, vm[i] == bus["vm"])

        #     for j in ref[:bus_gens][i]
        #         @constraint(model, pg[j] == ref[:gen][j]["pg"])
        #     end
        # end
    end

    optimize!(model)
    @info("LEM powerflow terminated with status: $(termination_status(model))")
    # Extract lineflow results
    p_ij = Dict()
    q_ij = Dict()

    for (l, pl) in p
        p_ij[l] = JuMP.value(pl)
    end
    for (l, ql) in q
        q_ij[l] = JuMP.value(ql)
    end

    return model, p_ij, q_ij, JuMP.value.(vm), JuMP.value.(va), JuMP.value.(pg), JuMP.value.(qg)
end
