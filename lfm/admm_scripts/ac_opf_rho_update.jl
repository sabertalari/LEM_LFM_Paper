function update_rho(ref, f_pl, f_ql, vm, va, f_plᵢ, f_qlᵢ, vmₗ, vaₗ, ρ_pl,ρ_ql,ρ_vm,ρ_va, μ_incr, μ_decr, τ_incr, τ_decr, ρ̲, ρ̅, k)

    for (i, bus) in ref[:bus] 
        for a in ref[:bus_arcs][i]

            r_pl = norm(f_plᵢ[string(a),k+1] - f_pl[string(a),k+1]) # pl[string(a),k+1] - plᵢ[string(a),k+1])
            s_pl = norm(f_pl[string(a),k+1] - f_pl[string(a),k]) # norm(plᵢ[string(a),k+1] - plᵢ[string(a),k])
            if r_pl > μ_incr * s_pl
                ρ_pl[string(a),k+1] = ρ_pl[string(a),k] * (1 + τ_incr)

            elseif s_pl > μ_decr * r_pl
                ρ_pl[string(a),k+1] = ρ_pl[string(a),k] * (1 + τ_decr)^-1

            else
                ρ_pl[string(a),k+1] = ρ_pl[string(a),k]
            end

            r_ql = norm(f_qlᵢ[string(a),k+1] - f_ql[string(a),k+1]) # ql[string(a),k+1] - qlᵢ[string(a),k+1])
            s_ql = norm(f_ql[string(a),k+1] - f_ql[string(a),k]) # norm(qlᵢ[string(a),k+1] - qlᵢ[string(a),k])
            if r_ql > μ_incr * s_ql
                ρ_ql[string(a),k+1] = ρ_ql[string(a),k] * (1 + τ_incr)

            elseif s_ql > μ_decr * r_ql
                ρ_ql[string(a),k+1] = ρ_ql[string(a),k] * (1 + τ_decr)^-1

            else
                ρ_ql[string(a),k+1] = ρ_ql[string(a),k]
            end

            r_vm = norm(vm[string(i),k+1] - vmₗ[string(a),k+1])
            s_vm = norm(vmₗ[string(a),k+1] - vmₗ[string(a),k])
            if r_vm > μ_incr * s_vm
                ρ_vm[string(a),k+1] = ρ_vm[string(a),k] * (1 + τ_incr)

            elseif s_vm > μ_decr * r_vm
                ρ_vm[string(a),k+1] = ρ_vm[string(a),k] * (1 + τ_decr)^-1

            else
                ρ_vm[string(a),k+1] = ρ_vm[string(a),k]
            end

            r_va = norm(va[string(i),k+1] - vaₗ[string(a),k+1])
            s_va = norm(vaₗ[string(a),k+1] - vaₗ[string(a),k])
            if r_va > μ_incr * s_va
                ρ_va[string(a),k+1] = ρ_va[string(a),k] * (1 + τ_incr)

            elseif s_va > μ_decr * r_va
                ρ_va[string(a),k+1] = ρ_va[string(a),k] * (1 + τ_decr)^-1

            else
                ρ_va[string(a),k+1] = ρ_va[string(a),k]
            end

        end
    end
    
    for i in keys(ρ_pl[:,k+1])
        if ρ_pl[i,k+1] > ρ̅
            ρ_pl[i,k+1] = ρ̅
        elseif ρ_pl[i,k+1] < ρ̲
            ρ_pl[i,k+1] = ρ̲
        end
    end
    for i in keys(ρ_ql[:,k+1])
        if ρ_ql[i,k+1] > ρ̅
            ρ_ql[i,k+1] = ρ̅
        elseif ρ_ql[i,k+1] < ρ̲
            ρ_ql[i,k+1] = ρ̲
        end
    end
    for i in keys(ρ_vm[:,k+1])
        if ρ_vm[i,k+1] > ρ̅
            ρ_vm[i,k+1] = ρ̅
        elseif ρ_vm[i,k+1] < ρ̲
            ρ_vm[i,k+1] = ρ̲
        end
    end
    for i in keys(ρ_va[:,k+1])
        if ρ_va[i,k+1] > ρ̅
            ρ_va[i,k+1] = ρ̅
        elseif ρ_va[i,k+1] < ρ̲
            ρ_va[i,k+1] = ρ̲
        end
    end
    return ρ_pl[:,k+1],ρ_ql[:,k+1],ρ_vm[:,k+1],ρ_va[:,k+1]
end