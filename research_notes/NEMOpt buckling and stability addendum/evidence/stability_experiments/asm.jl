using LinearAlgebra, SparseArrays
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

