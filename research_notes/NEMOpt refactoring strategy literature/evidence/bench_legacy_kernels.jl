# Evidence script for the NEMOpt refactoring study (research notes, not project code).
#
# Read-only probe of the legacy Full_Domain implementation: it parses the
# UNMODIFIED Full_Domain/NEMOpt.jl, evaluates the user's own function
# definitions, wraps them with call/factorization counters, and measures them.
# The numbers cited as "coordinator measurement" in reproduction_baseline.md
# came from this script.
#
# Environment used (2026-10-12): Julia 1.11.7 (linux x86_64, 4-core container),
# packages Tullio and Tensors (stdlib: SparseArrays, LinearAlgebra, Printf, Random).
# Run from the repository root:
#   NEL=40 julia --project=<env with Tullio, Tensors> "research_notes/NEMOpt refactoring strategy literature/evidence/bench_legacy_kernels.jl"
#   NEL=80 FULL=0 julia ...   # assembly/factorization scaling only
# Read-only probe of NEMOpt Full_Domain: evaluates the user's own function
# definitions (parsed from the unmodified source) and measures them.
using SparseArrays, LinearAlgebra, Tullio, Printf, Random

const SRC = normpath(joinpath(@__DIR__, "..", "..", "..", "Full_Domain", "NEMOpt.jl"))
const NEL = parse(Int, get(ENV, "NEL", "40"))

defname(x) = x isa Symbol ? x : (x isa Expr && x.head in (:call, :where, :(::))) ? defname(x.args[1]) : nothing
wanted = Set([:unitStrain, :convertVoigt, :H, :H_Deriv, :H_Steep, :H_Steep_Deriv, :g, :g_Deriv,
    :gaussPoints2D, :Generate_Mesh, :Filter_Shape_Functions, :constitutiveEqs, :Displacement_Solver,
    :χ_Solver, :Helmholtz_Filter, :Volume_Fraction, :Homogenizer_VolumeIntegral, :Macroscopic_Quantities,
    :Stiffness_Sensitivity_Filtered, :Stiffness_Sensitivity_Unfiltered, :Generate_Densities])
for ex in Meta.parseall(read(SRC, String)).args
    ex isa Expr || continue
    if ex.head in (:function, :(=)) && ex.args[1] isa Expr && ex.args[1].head in (:call, :where)
        if defname(ex.args[1]) == :constitutiveEqs
            ex.args[1].args[1] = :constitutiveEqs_orig   # load the original kernel under another name
        end
        (defname(ex.args[1]) in wanted || defname(ex.args[1]) == :constitutiveEqs_orig) && Core.eval(Main, ex)
    elseif ex.head == :(=) && ex.args[1] in (:e,) || (ex.head == :(=) && ex.args[1] isa Expr && ex.args[1].head == :tuple && :mapI in ex.args[1].args)
        Core.eval(Main, ex)
    end
end

convertVoigt(i, j) = intersect(findall(x -> x == i, mapI), findall(x -> x == j, mapJ))[1]

# Parameters exactly as in NEMOpt_main.jl (except mesh size)
Lx, Ly = 2.0, 2.0
Nelx, Nely = NEL, NEL
E, v = 100, 0.3
F_M = [1.025 0.00; 0.00 1.025]
r = 1.0
p, ϵ = 1, 5e-4
β, ω = 1.0, 0.5
β₂, ω₂ = 200.0, 0.05
gps = gaussPoints2D(2)

N_Elements, Ns, NodeNumbers, NC, XY, elCoords, DOFs, DOFs_Filter, D1, D2, D1_Filter = Generate_Mesh(Lx, Ly, Nelx, Nely)
println("mesh $(Nelx)x$(Nely): elements=$N_Elements, NDOF=$(Ns[4])")

# ---------- counters (wrap, do not modify, the original kernels) ----------
const NCONST = Ref(0)
const NSOLVE = Ref(0)
const _ce = constitutiveEqs_orig
constitutiveEqs(U, c, ξ, η, ρ) = (NCONST[] += 1; constitutiveEqs_orig(U, c, ξ, η, ρ))
Base.:\(A::SparseMatrixCSC{Float64,Int}, b::AbstractVecOrMat{Float64}) = (NSOLVE[] += 1; factorize(A) \ b)

# ---------- 1. single constitutive call ----------
Ue = zeros(8, 1); ce = elCoords[1, :, :]
constitutiveEqs(Ue, ce, 0.1, 0.2, 0.7)
t = @elapsed for _ in 1:2000; _ce(Ue, ce, 0.1, 0.2, 0.7); end
a = @allocated _ce(Ue, ce, 0.1, 0.2, 0.7)
@printf("constitutiveEqs: %.1f µs/call, %d bytes allocated/call\n", 1e6 * t / 2000, a)

# ---------- 2. filter undershoot ----------
Random.seed!(1)
ρbin = Float64.(rand(N_Elements) .> 0.5)
ρ̃, K_Filter = Helmholtz_Filter(ρbin, r * Lx / Nelx)
@printf("Helmholtz filter of random 0/1 design: min(nodal ρ̃) = %.3e, max = %.4f\n", minimum(ρ̃), maximum(ρ̃))
ρhalf = [elCoords[el, 1, 1] < 0 ? 1.0 : 0.0 for el in 1:N_Elements]
ρ̃h, _ = Helmholtz_Filter(ρhalf, r * Lx / Nelx)
ρone = zeros(N_Elements); ρone[div(N_Elements, 2) + div(Nelx, 2)] = 1.0
ρ̃o, _ = Helmholtz_Filter(ρone, r * Lx / Nelx)
@printf("  half solid/half void: min %.3e | single solid element: min %.3e\n", minimum(ρ̃h), minimum(ρ̃o))
println("  (-1e-4)^1.1 -> ", try (-1e-4)^1.1 catch err typeof(err) end)
println("Rep_MU for F_M=0.975I: ", abs(maximum([0.975 0; 0 0.975] - I)))

# ---------- 3. assembly: incremental CSC insertion vs triplets ----------
U0 = zeros(Float64, Ns[4], 1)
ρ̃s = clamp.(ρ̃, 0.0, 1.0)
function assemble_incremental(U)
    K = spzeros(Float64, Ns[4], Ns[4])
    for el = 1:N_Elements
        elK = spzeros(Float64, 8, 8)
        Ce = vec(DOFs[el, :, :]'); CeF = vec(DOFs_Filter[el, :, :]')
        for gp in eachrow(gps)
            ξ, η = gp[2:end]
            NF, _, _ = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
            ρx = (NF*ρ̃s[CeF])[1]
            B, detJ, _, _, _, A_rho_2D, _ = _ce(U[Ce, :], elCoords[el, :, :], ξ, η, ρx)
            elK += gp[1] * (B' * A_rho_2D * B) * detJ
        end
        K[Ce, Ce] += elK
    end
    K
end
function assemble_triplets(U)
    Is = Int[]; Js = Int[]; Vs = Float64[]
    sizehint!(Is, 64N_Elements); sizehint!(Js, 64N_Elements); sizehint!(Vs, 64N_Elements)
    for el = 1:N_Elements
        elK = zeros(8, 8)
        Ce = vec(DOFs[el, :, :]'); CeF = vec(DOFs_Filter[el, :, :]')
        for gp in eachrow(gps)
            ξ, η = gp[2:end]
            NF, _, _ = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
            ρx = (NF*ρ̃s[CeF])[1]
            B, detJ, _, _, _, A_rho_2D, _ = _ce(U[Ce, :], elCoords[el, :, :], ξ, η, ρx)
            elK .+= gp[1] * (B' * A_rho_2D * B) * detJ
        end
        for j in 1:8, i in 1:8
            push!(Is, Ce[i]); push!(Js, Ce[j]); push!(Vs, elK[i, j])
        end
    end
    sparse(Is, Js, Vs, Ns[4], Ns[4])
end
assemble_incremental(U0); assemble_triplets(U0)
t_inc = @elapsed K1 = assemble_incremental(U0)
t_tri = @elapsed K2 = assemble_triplets(U0)
@printf("tangent assembly: incremental K[Ce,Ce]+=elK %.2f s | triplets+sparse() %.2f s | same matrix: %s\n",
    t_inc, t_tri, isapprox(K1, K2))

# ---------- 4. which factorization does `\` pick for the reduced matrix? ----------
K = K2
Kff = K[1:end-2, 1:end-2]
m = 2 * Ns[1]; s = 2 * (Ns[1] + Ns[2]); n = 2 * (Ns[1] + Ns[2] + Ns[3] - 1)
Kmm, Kms, Kmi = Kff[1:m, 1:m], Kff[1:m, m+1:s], Kff[1:m, s+1:n]
Ksm, Kss, Ksi = Kff[m+1:s, 1:m], Kff[m+1:s, m+1:s], Kff[m+1:s, s+1:n]
Kim, Kis, Kii = Kff[s+1:n, 1:m], Kff[s+1:n, m+1:s], Kff[s+1:n, s+1:n]
Kred = vcat(hcat(Kmm + Kms * D1 + D1' * Ksm + D1' * Kss * D1, Kmi + D1' * Ksi), hcat(Kim + Kis * D1, Kii))
println("K symmetric: ", issymmetric(K), " | reduced K symmetric: ", issymmetric(Kred),
    " | asym norm: ", norm(Kred - Kred', Inf), " | factorize picks: ", nameof(typeof(factorize(Kred))))
rhs = randn(size(Kred, 1), 4)
factorize(Kred); cholesky(Symmetric(Kred))
t_lu = @elapsed for k in 1:4; factorize(Kred) \ rhs[:, k]; end
t_ch = @elapsed (Fc = cholesky(Symmetric(Kred)); Fc \ rhs)
@printf("4 RHS: 4x factorize+solve %.3f s | 1x cholesky + 4-RHS solve %.3f s\n", t_lu, t_ch)

# ---------- 5. one full objective+gradient evaluation, instrumented ----------
if get(ENV, "FULL", "1") == "1"
    ρd = rand(Int(N_Elements / 4))
    NCONST[] = 0; NSOLVE[] = 0
    t_obj = @elapsed begin
        ρu = Generate_Densities(ρd)
        global ρf, KF = Helmholtz_Filter(ρu, r * Lx / Nelx)
        n0 = NCONST[]; s0 = NSOLVE[]
        global Un = Displacement_Solver(ρf, 1e-7, 0.1)
        global n_newton = NCONST[] - n0; global s_newton = NSOLVE[] - s0
        global χ = [χ_Solver(mapI[i], mapJ[i], Un, ρf) for i in 1:4]
        n1 = NCONST[]
        P, A = Homogenizer_VolumeIntegral(Un, χ, ρf)
        global n_hom = NCONST[] - n1
    end
    ngp = 4 * N_Elements
    @printf("objective: %.1f s | constitutive calls/(elements*gauss pts): Newton %.0f, homogenization %.0f | factorizations: Newton %d, objective total %d\n",
        t_obj, n_newton / ngp, n_hom / ngp, s_newton, NSOLVE[])
    n2 = NCONST[]; s2 = NSOLVE[]
    t_grad = @elapsed Stiffness_Sensitivity_Unfiltered([(1, 1, 1, 1), (2, 2, 2, 2), (1, 1, 2, 2)], Un, χ, ρf, KF)
    @printf("gradient (3 stiffness sensitivities): %.1f s | constitutive calls/(el*gp) %.0f | factorizations %d\n",
        t_grad, (NCONST[] - n2) / ngp, NSOLVE[] - s2)
end

# ---------- 6. same constitutive law with Tensors.jl (type-stable, AD tangent) ----------
using Tensors
const λc, μc = E * v / ((1 + v) * (1 - 2v)), E / (2(1 + v))
ψ(F) = (J = det(F); μc / 2 * (tr(tdot(F)) - 2) - μc * log(J) + λc / 2 * log(J)^2)
stress_tangent(F) = ((∂²ψ, ∂ψ) = hessian(ψ, F, :all); (∂ψ, ∂²ψ))
Ur = 0.02 .* randn(8, 1)
B, _, P1, A2 = _ce(Ur, ce, 0.1, 0.2, 1.0)
NF, BF, _ = Filter_Shape_Functions(ce, 0.1, 0.2)
Fm = reshape(Ur, 2, 4) * BF' + I
Ft = Tensor{2,2}(Fm)
Pt, At = stress_tangent(Ft)
okP = all(isapprox(P1[q], Pt[mapI[q], mapJ[q]]; rtol=1e-10) for q in 1:4)
okA = all(isapprox(A2[a, b], At[mapI[a], mapJ[a], mapI[b], mapJ[b]]; rtol=1e-10, atol=1e-10) for a in 1:4, b in 1:4)
stress_tangent(Ft)
nrep = 10^6
function loop_st(Ft, n)
    acc = 0.0
    for _ in 1:n
        Pq, Aq = stress_tangent(Ft); acc += Pq[1, 1] + Aq[1, 1, 1, 1]
    end
    acc
end
loop_st(Ft, 10)
t = @elapsed loop_st(Ft, nrep)
@printf("Tensors.jl neo-Hookean P + tangent via AD: %.3f µs/call, %d bytes/call | matches original P: %s, tangent: %s\n",
    1e6 * t / nrep, @allocated(stress_tangent(Ft)), okP, okA)
