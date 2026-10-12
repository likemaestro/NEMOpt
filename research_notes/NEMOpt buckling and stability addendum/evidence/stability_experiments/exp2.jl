# Inertia reliability of CHOLMOD unpivoted LDL' vs exact eigenvalues, on a heterogeneous periodic cell.
using LinearAlgebra, SparseArrays, Printf, Random
include(joinpath(@__DIR__, "asm.jl"))
rng = MersenneTwister(7)
n = 14
dens = 0.02 .+ 0.98 .* rand(rng, n, n)          # strong contrast, no square symmetry
K = assemble(n; dens=dens)
ev = eigvals(Symmetric(Matrix(K)))
N = size(K,1)
nbad = 0; ntests = 0; nchol = 0
for m in vcat(1:10, 20:20:N-20)
    gap = ev[m+1] - ev[m]
    gap < 1e-8*abs(ev[m+1]) && continue
    σ = (ev[m] + ev[m+1])/2
    L = ldlt(Symmetric(K); shift=-σ, check=false)
    F = cholesky(Symmetric(K); shift=-σ, check=false)
    issuccess(F) && (nchol += 1)
    if issuccess(L)
        d = diag(sparse(L.LD)); neg = count(<(0), d)
        global ntests += 1
        neg != m && (global nbad += 1; @printf("  mismatch m=%d ldlt=%d\n", m, neg))
    else
        @printf("  ldlt failed at m=%d\n", m)
    end
end
@printf("heterogeneous n=%d N=%d: %d shifts tested, %d LDL' inertia mismatches, cholesky succeeded on %d indefinite shifts\n", n, N, ntests, nbad, nchol)
# near-coincident shift: sigma within 1e-9 relative of an eigenvalue
for rel in (1e-4, 1e-8, 1e-11, 1e-13)
    σ = ev[5]*(1 + rel)    # just above ev[5] -> exactly 5 negatives
    L = ldlt(Symmetric(K); shift=-σ, check=false)
    neg = issuccess(L) ? count(<(0), diag(sparse(L.LD))) : -1
    bk = LinearAlgebra.inertia(bunchkaufman(Symmetric(Matrix(K) - σ*I)))
    @printf("  sigma = ev5*(1+%g): ldlt #neg=%d (exact 5), bunchkaufman inertia=%s\n", rel, neg, string(bk))
end
# Cholesky failure position on an indefinite shift
σ = (ev[3] + ev[4])/2
try
    cholesky(Symmetric(K); shift=-σ)
catch e
    println("  cholesky throws ", typeof(e), ": ", sprint(showerror, e))
end
# Acoustic tensor of the FEAD2026 energy W = mu/2(I1-2) - mu lnJ + lam/2 (lnJ)^2, closed form vs scan
E, nu = 100.0, 0.3; lam = E*nu/((1+nu)*(1-2nu)); mu = E/(2(1+nu))
Atan(F) = begin
    Fi = inv(F)'; J = det(F); A = zeros(2,2,2,2)
    for i in 1:2, j in 1:2, k in 1:2, l in 1:2
        A[i,j,k,l] = mu*(i==k)*(j==l) + lam*Fi[i,j]*Fi[k,l] + (mu - lam*log(J))*Fi[i,l]*Fi[k,j]
    end; A end
function minQ(A; nθ=3600)
    best = Inf
    for θ in range(0, π, length=nθ+1)[1:end-1]
        nv = (cos(θ), sin(θ)); Q = zeros(2,2)
        for i in 1:2, j in 1:2, k in 1:2, l in 1:2; Q[i,k] += nv[j]*A[i,j,k,l]*nv[l]; end
        best = min(best, eigmin(Symmetric(Q)))
    end; best end
closed(F) = begin J = det(F); s = svdvals(F); c = lam + mu - lam*log(J)
    min(mu, mu + c/minimum(s)^2, mu + c/maximum(s)^2) end
for (lbl, F) in (("equibiaxial 0.5", [0.5 0; 0 0.5]), ("uniaxial 0.3", [0.3 0; 0 1.0]), ("equibiaxial 2.0", [2.0 0; 0 2.0]),
                 ("uniaxial 12", [12.0 0; 0 1.0]), ("shear 0.8", [1.0 0.8; 0 1.0]), ("biaxial 6x3", [6.0 0; 0 3.0]))
    @printf("  %-16s scan min eig Q = %+.6e  closed form = %+.6e\n", lbl, minQ(Atan(F)), closed(F))
end
# uniaxial tension stretch at loss of ellipticity: lam(ln s - 1) - mu = mu  (smallest stretch = 1)
s_crit = exp(1 + 2mu/lam); @printf("  uniaxial-tension LoE stretch (closed form) = %.4f; check minQ just above: %+.3e\n", s_crit, minQ(Atan([s_crit*1.001 0; 0 1.0])))
t = @elapsed for _ in 1:100; minQ(Atan([0.9 0.1; 0.05 0.95]); nθ=360); end
@printf("  cost of 360-angle scan (naive, allocating): %.2e s per state\n", t/100)
