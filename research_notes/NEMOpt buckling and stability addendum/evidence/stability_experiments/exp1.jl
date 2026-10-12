# Small API-behavior and cost experiment for the stability addendum (stdlib only).
using LinearAlgebra, SparseArrays, Printf, Random
BLAS.set_num_threads(1)
println("Julia ", VERSION, "  SuiteSparse/CHOLMOD ", SparseArrays.CHOLMOD.BUILD_VERSION)
println(methods(LinearAlgebra.inertia))

# ---------- periodic Q4 plane-strain linear-elastic stiffness, optional Bloch phases ----------
function q4ke(E, nu, h)
    lam = E*nu/((1+nu)*(1-2nu)); mu = E/(2(1+nu))
    C = [lam+2mu lam 0; lam lam+2mu 0; 0 0 mu]
    Ke = zeros(8,8); g = 1/sqrt(3)
    for (xi,eta) in ((-g,-g),(g,-g),(g,g),(-g,g))
        dNdxi = 0.25*[-(1-eta) (1-eta) (1+eta) -(1+eta)]
        dNdeta = 0.25*[-(1-xi) -(1+xi) (1+xi) (1-xi)]
        J = h/2; dNdx = dNdxi/J; dNdy = dNdeta/J
        B = zeros(3,8)
        for a in 1:4
            B[1,2a-1] = dNdx[a]; B[2,2a] = dNdy[a]; B[3,2a-1] = dNdy[a]; B[3,2a] = dNdx[a]
        end
        Ke .+= B'*C*B*J^2
    end
    Ke
end
# n x n periodic cell of side 1; unique nodes (i,j) i,j=0..n-1; returns K(k) with node `fix` removed (or none)
function assemble(n; k=(0.0,0.0), fixcenter=true, dens=nothing)
    h = 1.0/n; Ke = q4ke(100.0, 0.3, h)
    T = (k[1]==0 && k[2]==0) ? Float64 : ComplexF64
    I_ = Int[]; J_ = Int[]; V_ = T[]
    node(i,j) = mod(i,n) + n*mod(j,n) + 1
    for j in 0:n-1, i in 0:n-1
        corners = ((i,j),(i+1,j),(i+1,j+1),(i,j+1))
        dofs = Int[]; ph = T[]
        for (a,b) in corners
            nd = node(a,b); s = T(1)
            if T <: Complex
                s = cis(k[1]*(a>=n ? 1.0 : 0.0) + k[2]*(b>=n ? 1.0 : 0.0)) # phase e^{i k.L}, L=1
            end
            append!(dofs, (2nd-1, 2nd)); append!(ph, (s,s))
        end
        w = dens === nothing ? 1.0 : dens[i+1,j+1]
        for a in 1:8, b in 1:8
            push!(I_, dofs[a]); push!(J_, dofs[b]); push!(V_, w*conj(ph[a])*Ke[a,b]*ph[b])
        end
    end
    N = 2n^2
    K = sparse(I_, J_, V_, N, N)
    if fixcenter
        c = node(n÷2, n÷2); keep = setdiff(1:N, (2c-1, 2c))
        K = K[keep, keep]
    end
    K = (K + K')/2
    return K
end

# ---------- 1. API behavior on a small indefinite matrix ----------
n = 12
K = assemble(n)
ev = eigvals(Symmetric(Matrix(K)))
println("small n=$n  size=", size(K,1), "  lambda_min=", ev[1])
for m in (1, 3, 7)
    σ = (ev[m] + ev[m+1])/2   # exactly m eigenvalues below σ
    Ks = K - σ*I
    Ks = sparse(Symmetric(Ks))
    F = cholesky(Symmetric(Ks); check=false)
    L = ldlt(Symmetric(Ks); check=false)
    nneg = issuccess(L) ? count(<(0), diag(sparse(L.LD))) : -1
    Ls = ldlt(Symmetric(K); shift=-σ, check=false)
    nneg2 = issuccess(Ls) ? count(<(0), diag(sparse(Ls.LD))) : -1
    bk = bunchkaufman(Symmetric(Matrix(Ks)))
    inert = isdefined(LinearAlgebra, :inertia) ? LinearAlgebra.inertia(bk) : "n/a"
    @printf("m=%d: cholesky issuccess=%s  ldlt issuccess=%s  #neg(D)=%d  ldlt(shift) #neg=%d  bunchkaufman inertia=%s\n",
            m, issuccess(F), issuccess(L), nneg, nneg2, string(inert))
end
# exception type and info from a failing cholesky
try
    cholesky(Symmetric(sparse(Symmetric(K - ((ev[1]+ev[2])/2)*I))))
catch e
    println("cholesky(check=true) throws: ", typeof(e), " -> ", e)
end
# a matrix with a zero leading pivot in every symmetric ordering: [0 1; 1 0] (inertia (1,1,0))
A2 = sparse([0.0 1.0; 1.0 0.0])
L2 = ldlt(Symmetric(A2); check=false)
println("ldlt([0 1;1 0]) issuccess=", issuccess(L2))
# singular K without center fix (pure periodic, 2 rigid translations)
Kp = assemble(n; fixcenter=false)
Fp = cholesky(Symmetric(Kp); check=false); Lp = ldlt(Symmetric(Kp); check=false)
evp = eigvals(Symmetric(Matrix(Kp)))
println("no center fix: 3 smallest eig=", evp[1:3], "  cholesky issuccess=", issuccess(Fp), "  ldlt issuccess=", issuccess(Lp))
if issuccess(Lp); d = diag(sparse(Lp.LD)); println("   min |D| = ", minimum(abs, d), "  #neg=", count(<(0), d)); end
# k-zero limit with Bloch assembly at tiny k
Kk = assemble(n; k=(1e-6, 0.0), fixcenter=false)
evk = eigvals(Hermitian(Matrix(Kk)))
println("pure periodic, k=1e-6: 3 smallest eig=", evk[1:3])
Kk = assemble(n; k=(pi, 0.0), fixcenter=false)
println("pure periodic, k=(pi,0): 3 smallest eig=", eigvals(Hermitian(Matrix(Kk)))[1:3])
Kk = assemble(n; k=(pi, 0.0), fixcenter=true)
println("center fixed, k=(pi,0): 3 smallest eig=", eigvals(Hermitian(Matrix(Kk)))[1:3])

# ---------- 2. timing at 80x80 and 100x100 ----------
function lanczos_si(F, N, nev; maxit=60, T=Float64)
    # shift-invert Lanczos with full reorthogonalization on op x -> F \ x ; returns Ritz values of K^-1
    Q = zeros(T, N, maxit+1); α = zeros(maxit); β = zeros(maxit)
    q = randn(MersenneTwister(1), T, N); q ./= norm(q); Q[:,1] = q
    nsolve = 0; vals = Float64[]
    for j in 1:maxit
        w = F \ Q[:,j]; nsolve += 1
        α[j] = real(dot(Q[:,j], w))
        w .-= α[j] .* Q[:,j]; j > 1 && (w .-= β[j-1] .* Q[:,j-1])
        for _ in 1:2; w .-= Q[:,1:j]*(Q[:,1:j]'*w); end
        β[j] = norm(w); Q[:,j+1] = w ./ β[j]
        if j >= nev+5
            Tm = SymTridiagonal(α[1:j], β[1:j-1]); E = eigen(Tm)
            idx = sortperm(E.values, rev=true)[1:nev]
            res = abs.(β[j] .* E.vectors[end, idx])
            if maximum(res ./ abs.(E.values[idx])) < 1e-8
                return 1 ./ E.values[idx], nsolve
            end
        end
    end
    return vals, nsolve
end
tmin(f; r=3) = minimum(begin; t = @elapsed f(); t; end for _ in 1:r)
for n in (80, 100)
    K = assemble(n); N = size(K,1)
    S = Symmetric(K)
    tl = tmin(() -> lu(K))
    tc = tmin(() -> cholesky(S))
    F = cholesky(S)
    tnum = tmin(() -> cholesky!(F, S))
    tld = tmin(() -> ldlt(S))
    b = randn(N); ts = tmin(() -> F \ b; r=5)
    B4 = randn(N,4); ts4 = tmin(() -> F \ B4; r=5)
    nnzL = nnz(sparse(F.L))
    Kc = assemble(n; k=(0.7, 1.3)); Hc = Hermitian(Kc)
    tcc = tmin(() -> cholesky(Hc)); Fc = cholesky(Hc)
    tsc = tmin(() -> Fc \ complex.(b); r=5)
    tlz = @elapsed (vals, ns) = lanczos_si(F, N, 6)
    tlzc = @elapsed (valsc, nsc) = lanczos_si(Fc, N, 6; T=ComplexF64)
    @printf("n=%d N=%d nnz(K)=%d nnz(L)=%d | LU %.3f s | chol %.3f s | chol! numeric %.3f s | ldlt %.3f s | solve1 %.4f s | solve4 %.4f s\n",
            n, N, nnz(K), nnzL, tl, tc, tnum, tld, ts, ts4)
    @printf("      complex Hermitian chol %.3f s | complex solve %.4f s | SI-Lanczos real 6 eig: %d solves %.3f s | complex: %d solves %.3f s\n",
            tcc, tsc, ns, tlz, nsc, tlzc)
    println("      6 smallest eig (real, k=0, center fixed): ", round.(sort(vals); sigdigits=6))
    println("      6 smallest eig (complex k=(0.7,1.3)): ", round.(sort(valsc); sigdigits=6))
end
