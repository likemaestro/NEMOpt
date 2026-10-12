using LinearAlgebra, Random
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
