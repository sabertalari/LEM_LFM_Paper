function admm(ref, ϵ_abs, ϵ_rel, μ_decr, μ_incr, τ_incr, τ_decr, ρ̅, ρ̲, ρ_delay, fp_demand, fq_demand, p_pf, q_pf, vm_pf, va_pf, K)

    N_λ = 4 * length(ref[:arcs])
    # N_x = 2 * length(ref[:arcs]) + 2 * length(ref[:bus])
    # Index for named arrays
    arcs_index = [string(ref[:arcs][i]) for i = 1:length(ref[:arcs])]
    bus_index = [string(i) for i in keys(ref[:bus])]

    # Dual variables
    λ_pl = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(λ_pl, arcs_index, 1)
    λ_ql = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(λ_ql, arcs_index, 1)
    λ_vm = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(λ_vm, arcs_index, 1)
    λ_va = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(λ_va, arcs_index, 1)

    # Branch variables
    f_pl = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(f_pl, arcs_index, 1)
    f_ql = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(f_ql, arcs_index, 1)
    vmₗ = NamedArray(zeros(length(ref[:arcs]), K)) # Consensus
    setnames!(vmₗ, arcs_index, 1)
    vaₗ = NamedArray(zeros(length(ref[:arcs]), K)) # Consensus
    setnames!(vaₗ, arcs_index, 1)

    # Bus variables
    fp = NamedArray(zeros(length(ref[:gen]), K))
    fq = NamedArray(zeros(length(ref[:gen]), K))
    f_plᵢ = NamedArray(zeros(length(ref[:arcs]), K)) # Consensus
    setnames!(f_plᵢ, arcs_index, 1)
    f_qlᵢ = NamedArray(zeros(length(ref[:arcs]), K)) # Consensus
    setnames!(f_qlᵢ, arcs_index, 1)
    vm = NamedArray(zeros(length(ref[:bus]), K))
    setnames!(vm, bus_index, 1)
    va = NamedArray(zeros(length(ref[:bus]), K))
    setnames!(va, bus_index, 1)

    # cost = zeros(length(ref[:bus]), K)
    cost = NamedArray(zeros(length(ref[:bus]), K))
    setnames!(cost, bus_index, 1)
    dual_p_admm = zeros(length(ref[:bus]), K)
    dual_q_admm = zeros(length(ref[:bus]), K)

    for i in keys(ref[:bus])
        vm[string(i),1] = vm_pf[i]
        va[string(i),1] = va_pf[i]
    end

    for (l, i, j) in ref[:arcs]
        vmₗ[string((l, i, j)),1] = vm_pf[i]
        vaₗ[string((l, i, j)),1] = va_pf[i]


    end

    # Residuals
    r_norm = zeros(K)
    s_norm = zeros(K)

    # Primal residuals
    r_pl = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(r_pl, arcs_index, 1)
    r_ql = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(r_ql, arcs_index, 1)
    r_vm = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(r_vm, arcs_index, 1)
    r_va = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(r_va, arcs_index, 1)

    # Dual residuals
    s_pl = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(s_pl, arcs_index, 1)
    s_ql = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(s_ql, arcs_index, 1)
    s_vm = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(s_vm, arcs_index, 1)
    s_va = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(s_va, arcs_index, 1)

    # Feasibility tolerances for primal and dual residual
    ϵ_pri = zeros(K)
    ϵ_dual = zeros(K)

    # Penalty factors for each set of consensus variables
    ρ_pl = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(ρ_pl, arcs_index, 1)
    ρ_pl[:,1] .= ρ̲
    ρ_ql = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(ρ_ql, arcs_index, 1)
    ρ_ql[:,1] .= ρ̲
    ρ_vm = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(ρ_vm, arcs_index, 1)
    ρ_vm[:,1] .= ρ̲
    ρ_va = NamedArray(zeros(length(ref[:arcs]), K))
    setnames!(ρ_va, arcs_index, 1)
    ρ_va[:,1] .= ρ̲
    ρ = zeros(K)

    # Placeholder for last iteration, needed for plotting
    final_k = 0

    ### Start of ADMM loop ###
    for k in 1:K
        
        # Norm of penalty factor for evaluation during interations, not needed for ADMM
        ρ[k] = norm(cat(ρ_pl[:,k], ρ_ql[:,k], ρ_vm[:,k], ρ_va[:,k], dims=(1,1)))

        # Loop over all buses in the grid and individually update the buses
        for i in keys(ref[:bus])

            # a_lst = [string(a) for a in ref[:bus_arcs][i]]
            a_lst = []
            for a in ref[:bus_arcs][i]
                push!(a_lst,string(a))
            end
            
            fp[ref[:bus_gens][i],k+1], fq[ref[:bus_gens][i],k+1], f_plᵢ[a_lst,k+1], f_qlᵢ[a_lst,k+1], vm[string(i),k+1], va[string(i),k+1], dual_p_admm[i, k], dual_q_admm[i, k], cost[string(i), k] = update_single_bus(ref, vmₗ[:,k], vaₗ[:,k], f_pl[:,k], f_ql[:,k],
                                                                                                                                                                                                                        vm[string(i),k], va[string(i),k],
                                                                                                                                                                                                                        fp[:,k], fq[:,k], f_plᵢ[a_lst,k], f_qlᵢ[a_lst,k], 
                                                                                                                                                                                                                        P_lem, Q_lem, R_pᵈ, R_pᵘ, R_qᵈ, R_qᵘ,
                                                                                                                                                                                                                        λ_pl[:,k], λ_ql[:,k], λ_vm[:,k], λ_va[:,k], 
                                                                                                                                                                                                                        ρ_pl[:,k],ρ_ql[:,k],ρ_vm[:,k],ρ_va[:,k], i,
                                                                                                                                                                                                                        p_pf, q_pf,
                                                                                                                                                                                                                        fp_demand, fq_demand)
        end

        # Update all branches in the network
        for (l,i,j) in ref[:arcs]
            
            f_pl[[string((l,i,j)),string((l,j,i))],k+1], f_ql[[string((l,i,j)),string((l,j,i))],k+1], vmₗ[[string((l,i,j)),string((l,j,i))],k+1], vaₗ[[string((l,i,j)),string((l,j,i))],k+1] = update_single_branch(ref, vm[:,k+1], va[:,k+1], f_plᵢ[:,k+1], f_qlᵢ[:,k+1], 
                                                                                                                                                                                                                vmₗ[:,k], vaₗ[:,k], f_pl[:,k], f_ql[:,k], p_pf, q_pf, 
                                                                                                                                                                                                                λ_pl[:,k], λ_ql[:,k], λ_vm[:,k], λ_va[:,k], 
                                                                                                                                                                                                                ρ_pl[:,k],ρ_ql[:,k],ρ_vm[:,k],ρ_va[:,k],
                                                                                                                                                                                                                l, i, j)
        end

        # Update the dual variable, one for every pair of consensus variables
        λ_pl[:,k+1], λ_ql[:,k+1], λ_vm[:,k+1], λ_va[:,k+1] = update_lambda(
                                                                ref, λ_pl, λ_ql, λ_vm, λ_va, 
                                                                f_pl[:,k+1], f_ql[:,k+1], vm[:,k+1], va[:,k+1], 
                                                                f_plᵢ[:,k+1], f_qlᵢ[:,k+1], vmₗ[:,k+1], vaₗ[:,k+1], 
                                                                ρ_pl[:,k], ρ_ql[:,k], ρ_vm[:,k], ρ_va[:,k], k)

        # Update the primal and dual residuals for every pair of consensus variables, as well as their norm
        r_norm[k], s_norm[k], r_pl[:,k], r_ql[:,k], r_vm[:,k], r_va[:,k], s_pl[:,k], s_ql[:,k], s_vm[:,k], s_va[:,k]= update_residuals(ref, f_pl, f_ql, vm, va, f_plᵢ, f_qlᵢ, vmₗ, vaₗ, ρ_pl[:,k], ρ_ql[:,k], ρ_vm[:,k], ρ_va[:,k], r_pl, r_ql, r_vm, r_va, s_pl, s_ql, s_vm, s_va, k)  

        # Update the penalty factor, one for every pair of consensus variables
        if rem(k, ρ_delay) == 0
            ρ_pl[:,k+1],ρ_ql[:,k+1],ρ_vm[:,k+1],ρ_va[:,k+1] = update_rho(ref, f_pl, f_ql, vm, va, f_plᵢ, f_qlᵢ, vmₗ, vaₗ, ρ_pl, ρ_ql, ρ_vm, ρ_va, μ_incr, μ_decr, τ_incr, τ_decr, ρ̲, ρ̅, k)
        elseif k == 1
            ρ_pl[:,k+1],ρ_ql[:,k+1],ρ_vm[:,k+1],ρ_va[:,k+1] = update_rho(ref, f_pl, f_ql, vm, va, f_plᵢ, f_qlᵢ, vmₗ, vaₗ, ρ_pl, ρ_ql, ρ_vm, ρ_va, μ_incr, μ_decr, τ_incr, τ_decr, ρ̲, ρ̅, k)
        else 
            ρ_pl[:,k+1],ρ_ql[:,k+1],ρ_vm[:,k+1],ρ_va[:,k+1] = ρ_pl[:,k],ρ_ql[:,k],ρ_vm[:,k],ρ_va[:,k]
        end
        # Update the primal and dual feasibility tolerances
        ϵ_pri[k], ϵ_dual[k] = update_epsilons(ref, f_pl, f_ql, vmₗ, vaₗ, f_plᵢ, f_qlᵢ, vm, va, λ_pl, λ_ql, λ_vm, λ_va, N_λ, ϵ_abs, ϵ_rel, k)

        # Print result of current iteration every 20 iterations
        if rem(k, 20) == 0
            @info "Iter.: $(k) > ρ $(ρ[k]) r_norm $(r_norm[k]) ϵ_pri $(ϵ_pri[k])| s_norm $(s_norm[k]) ϵ_dual $(ϵ_dual[k]) | cost $(sum(cost[:,k]))"
        end

        # Check for convergence or end of max. iterations
        if r_norm[k] <= ϵ_pri[k] && s_norm[k] <= ϵ_dual[k]
            @info("ADMM terminates at iteration $(k)")
            final_k = k
            break
        elseif k == K-1  # -1 because otherwise updates above with k+1 run into an error in last iteration
            @info("Reached max iteration $(k)")
            final_k = k
            break
        elseif ϵ_dual[k] > 100
            @info("ADMM interrupted: ϵ_dual reached $(ϵ_dual[k])")
            final_k = k
            break
        end
    end

    return fp, fq, f_plᵢ, f_qlᵢ, vm, va, f_pl, f_ql, vmₗ, vaₗ, r_norm, s_norm, ϵ_pri, ϵ_dual, cost, dual_p_admm, dual_q_admm, final_k, ρ
        
end