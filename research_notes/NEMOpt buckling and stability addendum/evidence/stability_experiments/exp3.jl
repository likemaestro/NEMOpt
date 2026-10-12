# Boundary condensation for Bloch scans: inertia(K(θ)-σI) = inertia(Kii-σI) + inertia(S(θ)) (Haynsworth),
# with Kii k-independent. Check correctness on a small mesh; time per k-point at 80x80 and 100x100.
using LinearAlgebra, SparseArrays, Printf, Random
include(joinpath(@__DIR__, "asm.jl"))
BLAS.set_num_threads(1)
bnd(n) = sort(vcat([(2(i+1)-1, 2(i+1)) for i in 0:n-1]..., [(2(n*j+1)-1, 2(n*j+1)) for j in 1:n-1]...) |> x -> collect(Iterators.flatten(x)))
function split(n)
    b = bnd(n); N = 2n^2; inn = setdiff(1:N, b); return b, inn
end
inertia_dense(H) = (e = eigvals(Hermitian(Matrix(H))); (count(>(0), e), count(<(0), e)))
# ---- correctness on small mesh with heterogeneous density and a shift that makes some k unstable
n = 16; rng = MersenneTwister(3); dens = 0.05 .+ 0.95 .* rand(rng, n, n)
b, inn = split(n)
K0 = assemble(n; k=(1e-9, 0.0), fixcenter=false, dens=dens)   # complex path, ~k=0 (no fix)
Kii = real(K0[inn, inn])
lmin_ii = eigmin(Symmetric(Matrix(Kii)))
θs = range(-π, π, length=13)[1:end-1]
lmins = [eigmin(Hermitian(Matrix(assemble(n; k=(t1+1e-12, t2), fixcenter=false, dens=dens)))) for t1 in θs, t2 in θs]
σ = 0.5*(maximum(lmins) + minimum(lmins)); σ = min(σ, 0.9*lmin_ii)
@printf("small n=%d: lambda_min(Kii)=%.4g, lambda_min(K(θ)) in [%.4g, %.4g], shift σ=%.4g\n", n, lmin_ii, minimum(lmins), maximum(lmins), σ)
Fii = cholesky(Symmetric(Kii - σ*I))
nmis = 0; nunst = 0
for t1 in θs, t2 in θs
    K = assemble(n; k=(t1+1e-12, t2), fixcenter=false, dens=dens) - σ*I
    S = K[b,b] - K[b,inn]*(Fii \ Matrix(K[inn,b]))
    p1, m1 = inertia_dense(K); p2, m2 = inertia_dense(S)
    m1 > 0 && (global nunst += 1)
    m1 != m2 && (global nmis += 1)
end
@printf("  %d k-points, %d unstable (K(θ)-σI indefinite), %d inertia mismatches between full K and condensed S\n", length(θs)^2, nunst, nmis)
# ---- timing at 80 and 100: setup = factor Kii + solves for boundary couplings; per-k = dense S(θ) + Bunch-Kaufman/eigvals
for n in (80, 100)
    b, inn = split(n); nb = length(b)
    Kθ = assemble(n; k=(0.3, 0.2), fixcenter=false)
    Kii = real(Kθ[inn, inn])
    t_fac = @elapsed Fii = cholesky(Symmetric(Kii))
    # decompose K_ib(θ) into offset classes by comparing two phase settings (exact for this mesh): solve for each
    t_setup = @elapsed begin
        Kib = Kθ[inn, b]; Y = Fii \ Matrix(Kib)       # one θ only: nb complex RHS (the real cost scales x4 for 4 offset classes)
    end
    t_S = @elapsed S = Kθ[b,b] - Kθ[b,inn]*Y
    t_eig = @elapsed eigvals(Hermitian(Matrix(S)))
    t_bk = @elapsed bunchkaufman(Hermitian(Matrix(S)))
    t_cc = @elapsed cholesky(Hermitian(Kθ + 1e-9I))
    @printf("n=%d: N=%d boundary DOFs=%d | Kii chol %.3f s | %d-RHS solve %.3f s | dense S form %.3f s | dense eig(S) %.4f s | BK(S) %.4f s | full complex chol %.3f s\n",
            n, size(Kθ,1), nb, t_fac, nb, t_setup, t_S, t_eig, t_bk, t_cc)
end
# ---- clean shift-invert Lanczos timing at 80x80 (warm)
include(joinpath(@__DIR__, "lanczos.jl"))
for n in (80, 100)
    K = assemble(n); F = cholesky(Symmetric(K)); N = size(K,1)
    lanczos_si(F, N, 6)
    t = @elapsed (v, ns) = lanczos_si(F, N, 6)
    Kc = assemble(n; k=(0.7,1.3)); Fc = cholesky(Hermitian(Kc)); lanczos_si(Fc, N, 6; T=ComplexF64)
    tc = @elapsed (vc, nsc) = lanczos_si(Fc, N, 6; T=ComplexF64)
    @printf("n=%d warm SI-Lanczos (naive, full reorth.) real: %d solves %.3f s; complex: %d solves %.3f s\n", n, ns, t, nsc, tc)
end
