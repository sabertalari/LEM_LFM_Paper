function update_epsilons(ref, f_pl, f_ql, vmₗ, vaₗ, f_plᵢ, f_qlᵢ, vm, va, λ_pl, λ_ql, λ_vm, λ_va, N_λ, ϵ_abs, ϵ_rel, k)

    e_p = 0
    e_d = 0
    x = 0
    z = 0
    λ = 0

    x = norm(cat(f_plᵢ[:,k], f_qlᵢ[:,k], vm[:,k], va[:,k], dims=(1,1)))
    z = norm(cat(f_pl[:,k], f_ql[:,k], vmₗ[:,k], vaₗ[:,k], dims=(1,1)))
    λ = norm(cat(λ_pl[:,k], λ_ql[:,k], λ_vm[:,k], λ_va[:,k], dims=(1,1)))

    e_p = sqrt(N_λ) * ϵ_abs + ϵ_rel * max(x,z,0)
    e_d = sqrt(N_λ) * ϵ_abs + ϵ_rel * λ

    return e_p, e_d
end