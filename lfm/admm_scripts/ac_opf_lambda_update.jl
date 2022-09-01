function update_lambda(ref, λ_pl, λ_ql, λ_vm, λ_va, f_pl, f_ql, vm, va, f_plᵢ, f_qlᵢ, vmₗ, vaₗ, ρ_pl, ρ_ql, ρ_vm, ρ_va, k)

    for (i, bus) in ref[:bus] 
        for a in ref[:bus_arcs][i]
            λ_pl[string(a),k+1] = λ_pl[string(a),k] + ρ_pl[string(a)] * (f_plᵢ[string(a)] - f_pl[string(a)]) # (pl[string(a)] - plᵢ[string(a)])
            λ_ql[string(a),k+1] = λ_ql[string(a),k] + ρ_ql[string(a)] * (f_qlᵢ[string(a)] - f_ql[string(a)]) # (ql[string(a)] - qlᵢ[string(a)])
            λ_vm[string(a),k+1] = λ_vm[string(a),k] + ρ_vm[string(a)] * (vm[string(i)] - vmₗ[string(a)])
            λ_va[string(a),k+1] = λ_va[string(a),k] + ρ_va[string(a)] * (va[string(i)] - vaₗ[string(a)])

        end
    end
    return λ_pl[:,k+1], λ_ql[:,k+1], λ_vm[:,k+1], λ_va[:,k+1]
end