include("NEMOpt-DR.jl");

# NEMOpt-DR_main.jl
# v1.0.0
# author: Murat GÜVEN

ELAPSED_TIME = @elapsed begin
    Lx, Ly = 2.0, 2.0                                   # Horizontal and vertical dimensions of the microstructure.
    Nelx, Nely = 80, 80                               # Number of elements along the horizontal and vertical dimensions of the microstructure.

    E, v = 100, 0.3                                     # Modulus of elasticity and Poisson's ratio of the microstructure constituents.

    F_M = [1.025 0.00                                   # Applied macroscopic strain on the microstructure.
           0.00 1.025]

    r = 1.0                                             # Helmholtz filter radius in terms of element size (Lx/Nelx).
    p, ϵ = 1, 5.0e-4                                    # Penalty factor, Minimum density.
    β, ω = 1.0, 0.5                                     # Steepness of Heaviside function, Center of Heaviside function.
    β₂, ω₂ = 200.0, 0.05                                # Parameters for the Steep Heaviside function

    gps = gaussPoints2D(2)                              # Number of Gauss quadratures (n x n)

    N_ITER = 300                                        # Specify the maximum number of iterations for the optimization.
    V_MAX = 0.4                                         # Specify the maximum volume fraction for the topology.

    #RANDOM_SEED = 32218                                   # Specify the seed for rand() function.
    INITIAL_DESIGN_NAME = "ITER_1.png"                     # Specify the name of the topology that you want to feed.
    Description = ["Optimizing for negative Poisson ratio"]         # Describe the current analysis.
    global Jhist, Vfhist 
    Jhist = []
    Vfhist = []
    #---------------------------------------BEGINING OF THE ANALYSIS---------------------------------------
    println("Creating workspace... ")
    createWorkspace()

    println("Generating grid...")
    N_Elements, Ns, NodeNumbers, NC, XY, elCoords, DOFs, DOFs_Filter, D1, D2 = generateMesh(Lx, Ly, Nelx, Nely)
    NodeNumbers_Filter, Ns_Filter, DOFs_Filter, D1_Filter = generateMesh_Filter(Lx, Ly, Nelx, Nely)

    println("Creating Initial Design...")
    Initial_Design = Create_Initial_Design()

    function Objective(𝜌_Design)
        global 𝜌_Unfiltered, U_Nonlin, χ, 𝜌_Filtered, K_Filter, iter, Prev_J, p, β, P, A, v_12, v_0, Vf
        iter += 1

        begin # Unit Cell Analysis
            # Create element densities
            𝜌_Unfiltered = Generate_Densities(𝜌_Design)
            #println(𝜌_Unfiltered)
            # Apply the Helmholtz filter
            𝜌_Filtered, K_Filter = Helmholtz_Filter(𝜌_Unfiltered, r * Lx / Nelx)
            # Find the displacement field
            U_Nonlin = Displacement_Solver(𝜌_Filtered, 1e-7, 0.1)
            # Find the χ field
            χ = [χ_Solver(mapI[i], mapJ[i], U_Nonlin, 𝜌_Filtered) for i in 1:4]
            # Homogenize the stress and stiffness
            P, A = Homogenizer_VolumeIntegral(U_Nonlin, χ, 𝜌_Filtered)
            # Current Volume Fraction
            Vf = Volume_Fraction(𝜌_Filtered)
        end

        # Calculate the Macroscopic Quantities Using Macroscopic Stiffness Tensor   
        v_12, K, G = Macroscopic_Quantities(A)

        # Target Poisson's ratio
        v_0 = -1.0

        # Objective Function
        J = 100 * (v_12 - v_0)^2
        #J = 1 / (A[3, 3] + A[4, 4])
        #J = -100.0*1 / 4 * (A[1,1] + A[2,2] + A[1, 2] + A[2, 1])
        
        Plot_Filtered_Topology(𝜌_Filtered, iter)

        p, β = round(p, digits=2), round(β, digits=2)
        Change = round((J - Prev_J) / abs(Prev_J) * 100, digits=3)
        J, Vf = round(J, digits=3), round(Vf, digits=3)
        v_12, K, G = round(v_12, digits=3), round(K, digits=3), round(G, digits=3)
        row = [iter J Vf v_12 K G]
        push!(df, row)
        Prev_J = J

        writedlm("Homogenized_Stiffnesses/A_ITER_$iter.txt", round.(A, digits=4))
        writedlm("Homogenized_Stresses/P_ITER_$iter.txt", round.(P, digits=4))
        push!(Jhist,J)
        push!(Vfhist,Vf)
        println("Iter: $iter | Change: $Change% | J = $J, Vf = $Vf, v_12 = $v_12, K = $K, G = $G | p: $p, β: $β")

        # Increase p by 0.2 for every 5 iterations until it reaches 3
        if iter % 5 == 0 && p < 3
            p += 0.1
            p = round(p, digits=2)
        end
        # When p reaches 3 increase β by 1 for every 10 iterations until it reaches 40
        if p == 3
            if iter % 5 == 0 && β < 32
                β += 1
                β = round(β, digits=2)
                options = MMAOptions(s_init=0.5 / (β + 1), s_incr=1.2, s_decr=0.7, maxiter=N_ITER, maxinner=40, verbose=false)   #maxiter, outer_maxiter, maxinner,  verbose
            end
        end

        return J
    end

    function Objective_Grad(𝜌_Design)
        # Calculate stiffness sensitivities
        δA_1111, δA_2222, δA_1122 = Stiffness_Sensitivity_Unfiltered([(1, 1, 1, 1), (2, 2, 2, 2), (1, 1, 2, 2)], U_Nonlin, χ, 𝜌_Filtered, K_Filter)
        # Calculate stiffness sensitivities
        #δA_1212, δA_2121 = Stiffness_Sensitivity_Unfiltered([(1, 2, 1, 2), (2, 1, 2, 1)], U_Nonlin, χ, 𝜌_Filtered, K_Filter)

        # Derivative of v_12
        ∇v_12 = 2 * δA_1122 / (A[1, 1] + A[2, 2]) - 2 * A[1, 2] * (δA_1111 + δA_2222) / (A[1, 1] + A[2, 2])^2

        # Derivative of J
        ∇J = 200 * (v_12 - v_0) * ∇v_12

        # Derivative of J
        #∇J = -1 / (A[3, 3] + A[4, 4])^2 * (δA_1212 + δA_2121)
        # Derivative of J
        #∇J = -100.0*1 / 4* (δA_1111 + δA_2222 + δA_1122 + δA_2211)
        return ∇J
    end

    function Volume_Constraint(𝜌_Design)
        # Create element densities
        𝜌_Unfiltered = Generate_Densities(𝜌_Design)
        # Apply the Helmholtz filter
        𝜌_Filtered, K_Filter = Helmholtz_Filter(𝜌_Unfiltered, r * Lx / Nelx)

        Vf = Volume_Fraction(𝜌_Filtered)

        g = Vf - V_MAX
        return g
    end

    function Volume_Constraint_Grad(𝜌_Design)
        Vf_F = zeros(Float64, Ns_Filter[4], 1)
        for el = 1:N_Elements
            Ce_Filter = vec(DOFs_Filter[el, :, :]')
            elF = zeros(Float64, 4, 1)
            for gp in eachrow(gps)
                ξ, η = gp[2:end]
                N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]
                elF += gp[1] * (H_Deriv(𝜌_x) * N_Filter') * detJ_Filter
            end
            Vf_F[Ce_Filter] += 2 / (Lx * Ly) * elF
        end

        begin
            Kmm = K_Filter[1:1*Ns_Filter[1], 1:1*Ns_Filter[1]]
            Kms = K_Filter[1:1*Ns_Filter[1], 1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2])]
            Kmi = K_Filter[1:1*Ns_Filter[1], 1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3])]

            Ksm = K_Filter[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2]), 1:1*Ns_Filter[1]]
            Kss = K_Filter[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2]), 1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2])]
            Ksi = K_Filter[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2]), 1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3])]

            Kim = K_Filter[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3]), 1:1*Ns_Filter[1]]
            Kis = K_Filter[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3]), 1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2])]
            Kii = K_Filter[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3]), 1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3])]

            A_Upper = hcat(Kmm + Kms * D1_Filter + D1_Filter' * Ksm + D1_Filter' * Kss * D1_Filter, Kmi + D1_Filter' * Ksi)
            A_Lower = hcat(Kim + Kis * D1_Filter, Kii)
            A_Reduced = vcat(A_Upper, A_Lower)

            Fm = Vf_F[1:1*Ns_Filter[1], :]
            Fs = Vf_F[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2]), :]
            Fi = Vf_F[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3]), :]

            Bm = (Fm + D1_Filter'Fs)
            Bi = Fi
            B_Reduced = vcat(Bm, Bi)

            w_Reduced = A_Reduced \ B_Reduced

            w_Masters = w_Reduced[1:1*Ns_Filter[1], :]
            w_Slaves = D1_Filter * w_Masters
            w_Internals = w_Reduced[1*Ns_Filter[1]+1:end, :]

            w_Filter = vcat(w_Masters, w_Slaves, w_Internals)
        end

        Vf_δθ_Full_Unfiltered = zeros(Float64, N_Elements, 1)
        for el = 1:N_Elements
            el_Vf_δθ = 0.0
            Ce_Filter = vec(DOFs_Filter[el, :, :]')
            for gp in eachrow(gps)
                ξ, η = gp[2:end]
                N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                el_Vf_δθ += gp[1] * (N_Filter*w_Filter[Ce_Filter])[1] * detJ_Filter
            end
            Vf_δθ_Full_Unfiltered[el] = el_Vf_δθ
        end

        Vf_δθ_Unfiltered = zeros(Int(N_Elements / 2), 1)
        elements_BL = reduce(hcat, [vcat(Int(i * Nelx + 1):Int(i * Nelx + Nelx / 2)) for i = 0:Nely/2-1])'
        for k = 1:Int(N_Elements / 2)
            el = elements_BL[k]
            i, j = Int(ceil(el / Nely)), (el - 1) % Nely + 1
            elV = el + Nelx - (2 * j - 1)
            Vf_δθ_Unfiltered[k] = Vf_δθ_Full_Unfiltered[el] + Vf_δθ_Full_Unfiltered[elV]
        end

        ∇g = Vf_δθ_Unfiltered

        return ∇g
    end

    function ChainRulesCore.rrule(::typeof(Objective), 𝜌_Design::AbstractVector)
        J = Objective(𝜌_Design)
        ∇J = Objective_Grad(𝜌_Design)
        J, Δ -> (NoTangent(), ∇J * Δ)
    end

    function ChainRulesCore.rrule(::typeof(Volume_Constraint), 𝜌_Design::AbstractVector)
        g = Volume_Constraint(𝜌_Design)
        ∇g = Volume_Constraint_Grad(𝜌_Design)
        g, Δ -> (NoTangent(), ∇g * Δ)
    end

    println("Creating model...")
    model = Model(Objective)

    println("Creating variables...")
    addvar!(model, zeros(Int(N_Elements / 2)), ones(Int(N_Elements / 2)))

    println("Adding constraints...")
    add_ineq_constraint!(model, 𝜌 -> Volume_Constraint(𝜌))

    println("Assigning MMA parameters...")
    alg = MMA02() # or MMA87()
    options = MMAOptions(tol=Nonconvex.Tolerance(kkt=1e-8, f=0.0),
        s_init=0.25, s_incr=1.2, s_decr=0.7, maxiter=N_ITER, maxinner=40, verbose=false)   #maxiter, outer_maxiter, maxinner,  verbose

    println("Starting optimization...")
    result = optimize(model, alg, Initial_Design, options=options, convcriteria=KKTCriteria())
    𝜌_Optimum = result.minimizer  # decision variables

    println("Optimization ended.")
end

writedlm("Elapsed_Time.txt", ELAPSED_TIME)
datfile = open("vf_hist.txt","w")
writedlm(datfile,Vfhist)  
close(datfile)
datfile = open("obj_hist.txt","w")
writedlm(datfile,Jhist)  
close(datfile)

println("Exporting results...")
Export_Results()
println("Done!")
