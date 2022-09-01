function update_residuals(ref, f_pl, f_ql, vm, va, f_plᵢ, f_qlᵢ, vmₗ, vaₗ, ρ_pl, ρ_ql, ρ_vm, ρ_va, r_pl, r_ql, r_vm, r_va, s_pl, s_ql, s_vm, s_va, k)

    r = 0.0
    s = 0.0
    for (i, bus) in ref[:bus] 
        for a in ref[:bus_arcs][i]

            # Lines := z, buses := x

            r_pl[string(a),k] = f_plᵢ[string(a),k+1] - f_pl[string(a),k+1] # pl[string(a),k+1] - plᵢ[string(a),k+1]
            r_ql[string(a),k] = f_qlᵢ[string(a),k+1] - f_ql[string(a),k+1] # ql[string(a),k+1] - qlᵢ[string(a),k+1]
            r_vm[string(a),k] = vm[string(i),k+1] - vmₗ[string(a),k+1]
            r_va[string(a),k] = va[string(i),k+1] - vaₗ[string(a),k+1]



            s_pl[string(a),k] = ρ_pl[string(a)] * (f_pl[string(a),k+1] - f_pl[string(a),k])
            s_ql[string(a),k] = ρ_ql[string(a)] * (f_ql[string(a),k+1] - f_ql[string(a),k])
            s_vm[string(a),k] = ρ_vm[string(a)] * (vmₗ[string(a),k+1] - vmₗ[string(a),k])
            s_va[string(a),k] = ρ_va[string(a)] * (vaₗ[string(a),k+1] - vaₗ[string(a),k])
            
        end
    end

    r = norm(cat(r_pl[:,k], r_ql[:,k], r_vm[:,k], r_va[:,k], dims=(1,1)))
    s = norm(cat(s_pl[:,k], s_ql[:,k], s_vm[:,k], s_va[:,k], dims=(1,1)))

    return r, s, r_pl[:,k], r_ql[:,k], r_vm[:,k], r_va[:,k], s_pl[:,k], s_ql[:,k], s_vm[:,k], s_va[:,k]

end