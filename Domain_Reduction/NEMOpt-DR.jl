# NEMOpt-DR.jl
# v1.0.0
# author: Murat GÜVEN

print("Importing libraries... ")
using Nonconvex, ChainRulesCore, Dates, Random, Distributions, SparseArrays, LinearAlgebra, Plots, Colors, Images, JLD2, DataFrames, Printf, DelimitedFiles, Tullio, FFMPEG, PlotlyJS;
Nonconvex.@load MMA
println("Imported.")

HOME_DIR = @__DIR__
FOLDER_NAME = "DR-" * Dates.format(now(), "yymmdd-HHMMSS")

Parameters = []
function createWorkspace()
    mkdir(FOLDER_NAME)
    cd(FOLDER_NAME)
    mkdir("Homogenized_Stiffnesses")
    mkdir("Homogenized_Stresses")
    mkdir("Topologies")
    mkdir("Graphs")

    Parameters = ["Length(X): $Lx", "Length(Y): $Ly", "Number of Elements(X): $Nelx", "Number of Elements(Y): $Nely",
        "Modulus of Elasticity: $E", "Poisson's Ratio: $v", "Helmholtz Filter Radius(el): $r", "Penalty Factor: $p",
        "Minimum Density: $ϵ", "Deformation Gradient: $F_M", "Max Number of Iterations: $N_ITER",
        "Volume Fraction: $V_MAX", "β: $β", "ω: $ω", "β₂: $β₂", "ω₂: $ω₂"]

    if !xor((@isdefined INITIAL_DESIGN_NAME), (@isdefined RANDOM_SEED))
        throw("Please provide 'INITIAL_DESIGN_NAME' OR 'RANDOM_SEED'")
    elseif @isdefined INITIAL_DESIGN_NAME
        push!(Parameters, "Design: $INITIAL_DESIGN_NAME")
    elseif @isdefined RANDOM_SEED
        push!(Parameters, "Seed: $RANDOM_SEED")
    end

    writedlm("Parameters.txt", Parameters)
    writedlm("Description.txt", Description)
end
Plots.default(fontfamily="Computer Modern", linewidth=2, framestyle=:box)

mapI, mapJ = [1, 2, 1, 2], [1, 2, 2, 1]

e = [[1.0, 0.0] [0.0, 1.0]]
function unitStrain(k, l)
    return e[k, :] * e[l, :]'
end

convertVoigt(i, j) = intersect(findall(x -> x == i, mapI), findall(x -> x == j, mapJ))[1];

H(𝜌) = (tanh(β * ω) + tanh(β * (𝜌 - ω))) / (tanh(β * ω) + tanh(β * (1 - ω)))
H_Deriv(𝜌) = β * sech(β * (𝜌 - ω))^2 / (tanh(β * ω) + tanh(β * (1 - ω)))

H_Steep(𝜌) = (tanh(β₂ * ω₂) + tanh(β₂ * (𝜌 - ω₂))) / (tanh(β₂ * ω₂) + tanh(β₂ * (1 - ω₂)))
H_Steep_Deriv(𝜌) = β₂ * sech(β₂ * (𝜌 - ω₂))^2 / (tanh(β₂ * ω₂) + tanh(β₂ * (1 - ω₂)))

g(𝜌) = max(H(𝜌)^p, ϵ)
g_Deriv(𝜌) = H(𝜌)^p > ϵ ? p * H(𝜌)^(p - 1) * H_Deriv(𝜌) : 0.0

iter, Prev_J = 0, Inf
df = DataFrame(Iteration=Int64[], J=Float64[], Vf=Float64[], v_12=Float64[], K=Float64[], G=Float64[])

function Macroscopic_Quantities(A)
    # Poisson's ratio definition
    v_12 = 2 * A[1, 2] / (A[1, 1] + A[2, 2])
    # Bulk modulus definition
    K = 1 / 4 * (A[1, 1] + A[2, 2] + A[1, 2] + A[2, 1])
    # Shear modulus definition
    G = A[3, 3] 
    return v_12, K, G
end

function Volume_Fraction(𝜌_Filtered::Vector)
    Vf = 0.0
    for el = 1:N_Elements
        Ce_Filter = vec(DOFs_Filter[el, :, :]')
        Vf_x = 0.0
        for gp in eachrow(gps)
            ξ, η = gp[2:end]
            N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
            𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]
            Vf_x += gp[1] * H(𝜌_x) * detJ_Filter
        end
        Vf += 2 / (Lx * Ly) * Vf_x
    end
    return Vf
end

function Create_Initial_Design()
    # Create the initial design
    if !xor((@isdefined INITIAL_DESIGN_NAME), (@isdefined RANDOM_SEED))
        throw("Please provide 'INITIAL_DESIGN_NAME' OR 'RANDOM_SEED'")
    elseif @isdefined INITIAL_DESIGN_NAME
        initial_design_path = joinpath(HOME_DIR, "Initial_Designs", INITIAL_DESIGN_NAME)
        Initial_Design = getTopology(initial_design_path, true)
    elseif @isdefined RANDOM_SEED
        Random.seed!(RANDOM_SEED)
        Initial_Design = rand(Uniform(0.0, 1.0), Int(N_Elements / 2))
    end
    return Initial_Design
end

function Plot_Displacement(U) # FIX
    gr(size=(300, 300))

    xs = range(-Lx / 2, Lx / 2, length=Nelx + 1)
    ys = range(-Ly / 2, Ly / 2, length=Nely + 1)

    Part_L = U[NodeNumbers]
    Part_TL = Part_L[1:Int(Nely / 2 + 1), :]
    Part_BL = Part_L[Int(Nely / 2 + 2):end, :]

    Part_TR, Part_BR = Part_BL[:, 1:end-1], Part_TL[:, 1:end-1] #reverse(-Part_BL[:, 1:end-1], dims=2), reverse(-Part_TL[:,1:end-1], dims=2)
    Part_R = vcat(Part_TR, Part_BR)
    U_Full = hcat(Part_L, Part_R)

    global U_Full
    name = lpad(iter, 10, "0")
    figure = Plots.contourf(xs, ys, reverse(U_Full, dims=1), lw=0, levels=100, c=:thermal, colorbar=false, axis=nothing, showaxis=false, aspect_ratio=:equal)
#    Plots.savefig(figure, "Topologies/$name.png")
end

function Plot_Filtered_Topology(𝜌_Filtered, iter)
    gr(size=(500, 500))

    xs = range(-Lx / 2, Lx / 2, length=Nelx + 1)
    ys = range(-Ly / 2, Ly / 2, length=Nely + 1)
    ysII = range(-Ly / 2, 0.0, length=Int.(Nely/2)+1)

    Part_B = 𝜌_Filtered[NodeNumbers_Filter]
    Part_T = reverse(Part_B[2:end, 1:end], dims=1)
    𝜌_Filtered_2D = vcat(Part_T, Part_B)

    name = lpad(iter, 10, "0")
    nameII = lpad(iter, 10, "h")
    figure = Plots.contourf(xs, ys, H.(𝜌_Filtered_2D), lw=0, levels=50, clims=(0, 1), colorbar=false, axis=nothing, showaxis=false, c=cgrad([:white, :black]), aspect_ratio=:equal)
    figureII = Plots.contourf(xs, ysII, H.(Part_B), lw=0, levels=50, clims=(0, 1), colorbar=false, axis=nothing, showaxis=false, c=cgrad([:white, :black]), aspect_ratio=:equal)
    Plots.savefig(figure, "Topologies/$name.png")
    Plots.savefig(figure, "Topologies/$name.svg")
    Plots.savefig(figureII, "Topologies/$nameII.png")
    Plots.savefig(figureII, "Topologies/$nameII.svg")
end

function Plot_Tiled_Filtered_Topology(𝜌_Filtered, tile_size)
    gr(size=(500 * tile_size, 500 * tile_size))
    xs = range(-Lx / 2 * tile_size, Lx / 2 * tile_size, length=(Nelx) * tile_size)
    ys = range(-Ly / 2 * tile_size, Ly / 2 * tile_size, length=(Nelx) * tile_size)

    Part_B = 𝜌_Filtered[NodeNumbers_Filter]
    Part_T = reverse(Part_B[2:end, 1:end], dims=1)
    𝜌_Filtered_2D = H.(vcat(Part_T, Part_B))

    Tiled_Topology = repeat(𝜌_Filtered_2D[2:end, 1:end-1], outer=[tile_size, tile_size])

    figure = Plots.contourf(xs, ys, Tiled_Topology, lw=0, levels=50, clims=(0, 1), colorbar=false, axis=nothing, showaxis=false, c=cgrad([:white, :black]), aspect_ratio=:equal)
    Plots.savefig(figure, "TILED_($(tile_size)x$(tile_size)).png")
    Plots.savefig(figure, "TILED_($(tile_size)x$(tile_size)).svg")
end

function getTopology(path, normalize)
    img = load(path)
    img_Size = size(img)
    img_Cropped = img[floor(Int, img_Size[1] / 2):img_Size[1], 1:floor(Int, img_Size[1] / 2)]
    img_Resized = imresize(img_Cropped, floor(Int, Nely / 2), floor(Int, Nely / 2))
    img_Gray = channelview(Gray.(img_Resized))
    Data = vec(reverse(float.(1 .- img_Gray), dims=1))

    𝜌_Filtered, K_Filter = Helmholtz_Filter(Generate_Densities(Data), r * Lx / Nelx)
    Vf = Volume_Fraction(𝜌_Filtered)

    if normalize
        Normalized_Data = Data * min(V_MAX / Vf, 1)
    else
        Normalized_Data = Data
    end

    return Normalized_Data
end

function Generate_Densities(𝜌_Design)
    Quadrant_TL = reshape(𝜌_Design, Int(Nelx / 2), Int(Nely / 2))
    Quadrant_TR = reverse(Quadrant_TL, dims=2)
    Part_T = hcat(Quadrant_TL, Quadrant_TR)
    𝜌_Unfiltered = vec(Part_T')
    return 𝜌_Unfiltered
end

function Export_Results()
    Titles = ["Objective", "Volume Fraction", "Poisson's Ratio", "Bulk Modulus", "Shear Modulus"]
    gr(size=(800, 500))
    for i = 2:ncol(df)
        title = Titles[i-1]
        name = names(df)[i]
        figure = Plots.plot(Matrix(df)[:, i], labels="$name", legend=:topright, xlabel="Iteration", ylabel=title)
        Plots.savefig(figure, "Graphs/$title.png")
        Plots.savefig(figure, "Graphs/$title.svg")
    end
    𝜌_Unfiltered = Generate_Densities(𝜌_Optimum)
    𝜌_Filtered, K_Filter = Helmholtz_Filter(𝜌_Unfiltered, r * Lx / Nelx)

    Plot_Tiled_Filtered_Topology(𝜌_Filtered, 3)
    Plot_Tiled_Filtered_Topology(𝜌_Filtered, 5)

    framerate = 30 #min(ceil(N_ITER/5), 30)
    FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -i Topologies/%10d.png -vf "scale=500:500" -c:v libx264 -pix_fmt yuv420p -y "Optimum-Topology.mp4" -hide_banner -loglevel error`)

    #save_object("Graphs/Final_Topology.jld2", 𝜌_Unfiltered);
    #save_object("Graphs/Optimization_His  tory.jld2", df);
end

function gaussPoints2D(n::Int64)
    # Gauss points
    coords, weights, x, w = zeros(n, n, 2), zeros(n, n), [], []
    if n == 2
        x = [-1 / sqrt(3), 1 / sqrt(3)]
        w = [1.0, 1.0]
    elseif n == 3
        x = [-sqrt(3 / 5), 0, sqrt(3 / 5)]
        w = [5 / 9, 8 / 9, 5 / 9]
    elseif n == 4
        x = [-sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5)), -sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5)), sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5)), sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5))]
        w = [(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36]
    end
    for i in 1:n
        for j = 1:n
            coords[i, j, :] = [x[i], x[j]]
            weights[i, j] = w[i] * w[j]
        end
    end
    gpCoord = reshape(coords, n * n, 2)
    gpW = reshape(weights, n * n, 1)
    return hcat(gpW, gpCoord)
end

function generateMesh(Lx::Float64, Ly::Float64, Nelx::Int64, Nely::Int64)
    dx, dy = (Lx / Nelx, Ly / Nely)
    N_Elements = Int(Nelx * Nely / 2)

    Nx, Ny = (Nelx + 1, Int(Nely / 2) + 1)

    N_Nodes = Nx * Ny
    N_Edges = 2 * Int(Nx - 3) + 2 * Int(Ny - 2)

    N_Prescribed = 6
    N_Masters = Int(Nx - 3) + Int(Ny - 2)
    N_Slaves = Int(Nx - 3) + Int(Ny - 2)
    N_Internals = N_Nodes - (N_Masters + N_Slaves + N_Prescribed)

    NDOF = 2 * N_Nodes

    MidCol = ceil(Int, Nx / 2)

    NodeNumbers = zeros(Int64, Ny, Nx)
    # PRESCRIBED NODES
    NodeNumbers[1, 1], NodeNumbers[end, 1], NodeNumbers[1, MidCol], NodeNumbers[end, MidCol], NodeNumbers[1, end], NodeNumbers[end, end] = 1:N_Prescribed
    # MASTER NODES
    NodeNumbers[1, 2:Int(MidCol - 1)] = N_Prescribed+1:N_Prescribed+Int((Nx - 3) / 2) # EDGE TOP LEFT
    NodeNumbers[end, 2:Int(MidCol - 1)] = N_Prescribed+Int((Nx - 3) / 2)+1:N_Prescribed+Int(Nx - 3) # EDGE BOTTOM LEFT
    NodeNumbers[2:Int(Ny - 1), 1] = N_Prescribed+Int(Nx - 3)+1:N_Prescribed+N_Masters # EDGE LEFT
    # SLAVE NODES
    NodeNumbers[1, Int(MidCol + 1):Int(Nx - 1)] = N_Prescribed+N_Masters+Int((Nx - 3) / 2):-1:N_Prescribed+N_Masters+1 # EDGE TOP RIGHT
    NodeNumbers[end, Int(MidCol + 1):Int(Nx - 1)] = N_Prescribed+N_Masters+Int(Nx - 3):-1:N_Prescribed+N_Masters+Int((Nx - 3) / 2)+1 # EDGE BOTTOM RIGHT
    NodeNumbers[2:Int(Ny - 1), end] = N_Prescribed+N_Masters+Int(Nx - 3)+1:N_Prescribed+N_Edges # EDGE RIGHT
    # INTERNAL NODES
    NodeNumbers[2:end-1, 2:end-1] = N_Prescribed+N_Edges+1:N_Nodes

    # Generate connectivity of each element based on their node numbers.
    NC = zeros(Int64, N_Elements, 4)
    for j = 1:Int(Nely / 2)
        for i = 1:Nelx
            NC[i+Nelx*(j-1), :] = [NodeNumbers[end-j+1, i], NodeNumbers[end-j+1, i+1], NodeNumbers[end-j, i+1], NodeNumbers[end-j, i]]
        end
    end

    # Nodal Coordinates
    x, y = (-Lx/2:dx:Lx/2, 0:dy:Ly/2)
    XY = zeros(Float64, N_Nodes, 2)
    for j = 1:Ny
        for i = 1:Nx
            XY[i+Nx*(j-1), :] = [x[i], y[j]]
        end
    end

    # Element Coordinates
    elCoords = zeros(Float64, N_Elements, 4, 2)
    for j = 1:Int(Nely / 2)
        for i = 1:Nelx
            n = i + Nx * (j - 1)
            ind = [n, n + 1, n + 1 + Nx, n + Nx]
            elCoords[i+Nelx*(j-1), :, :] = XY[ind, :]
        end
    end

    # Degree of Freedoms of the Elements
    DOFs = reshape(hcat(2 * NC .- 1, 2 * NC), (N_Elements, 4, 2))

    # Dependency Matrix
    D1_rowInd = vcat(1:2*N_Slaves)
    D1_colInd = vcat(1:2*N_Masters)
    D1_data = vcat(-ones(2 * (Nx - 3)), ones(2 * (Ny - 2)))
    D1 = sparse(D1_rowInd, D1_colInd, D1_data, 2 * N_Slaves, 2 * N_Masters)

    D2_rowInd = vcat(1:(Nx-3), 2*(N_Slaves-(Ny-2))+1:2*N_Slaves)
    D2_colInd = vcat(repeat([3, 4], Int((Nx - 3) / 2)), repeat([1, 2], (Ny - 2)))
    D2_data = ones(Float64, 2 * (Ny - 2) + (Nx - 3))
    D2 = sparse(D2_rowInd, D2_colInd, D2_data, 2 * N_Slaves, 2 * 2)

    # Number of Degrees of Freedom for Helmholtz Filter
    NDOF_Filter = 1 * N_Nodes

    # Degree of Freedoms of the Elements for Helmholtz Filter
    DOFs_Filter = reshape(1 * NC, (N_Elements, 4, 1))

    Ns = [N_Prescribed, N_Masters, N_Slaves, N_Internals, NDOF, NDOF_Filter]

    println("Number of elements: $(N_Elements) ($Nelx x $Nely)")
    return N_Elements, Ns, NodeNumbers, NC, XY, elCoords, DOFs, DOFs_Filter, D1, D2
end

function generateMesh_Filter(Lx::Float64, Ly::Float64, Nelx::Int64, Nely::Int64)
    dx, dy = (Lx / Nelx, Ly / Nely)
    N_Elements = Int(Nelx * Nely / 2)

    Nx, Ny = (Nelx + 1, Int(Nely / 2) + 1)

    N_Nodes = Nx * Ny
    N_Edges = 2 * Int(Nx - 3) + 2 * Int(Ny - 2)

    N_Masters = Int(Nx - 3) + Int(Ny - 2) + 2
    N_Slaves = Int(Nx - 3) + Int(Ny - 2) + 2
    N_Internals = N_Nodes - (N_Masters + N_Slaves)

    NDOF = N_Nodes

    MidCol = ceil(Int, Nx / 2)

    NodeNumbers = zeros(Int64, Ny, Nx)
    # MASTER NODES
    NodeNumbers[1, 1] = 1
    NodeNumbers[end, 1] = 2
    NodeNumbers[1, 2:Int(MidCol - 1)] = 2+1:2+Int((Nx - 3) / 2) # EDGE TOP LEFT
    NodeNumbers[end, 2:Int(MidCol - 1)] = 2+Int((Nx - 3) / 2)+1:2+Int(Nx - 3) # EDGE BOTTOM LEFT
    NodeNumbers[2:Int(Ny - 1), 1] = 2+Int(Nx - 3)+1:N_Masters # EDGE LEFT
    # SLAVE NODES
    NodeNumbers[1, end] = N_Masters + 1
    NodeNumbers[end, end] = N_Masters + 2
    NodeNumbers[1, Int(MidCol + 1):Int(Nx - 1)] = N_Masters+2+Int((Nx - 3) / 2):-1:N_Masters+3 # EDGE TOP RIGHT
    NodeNumbers[end, Int(MidCol + 1):Int(Nx - 1)] = N_Masters+2+Int(Nx - 3):-1:N_Masters+2+Int((Nx - 3) / 2)+1 # EDGE BOTTOM RIGHT
    NodeNumbers[2:Int(Ny - 1), end] = N_Masters+2+Int(Nx - 3)+1:N_Masters+2+Int(Nx - 3)+Int(Ny - 2) # EDGE RIGHT
    # INTERNAL NODES
    NodeNumbers[2:end-1, 2:end-1] = N_Masters+2+Int(Nx - 3)+Int(Ny - 2)+1:N_Nodes-2
    NodeNumbers[1, MidCol] = N_Nodes - 1
    NodeNumbers[end, MidCol] = N_Nodes

    #NodeNumbers = reverse(NodeNumbers, dims=1)

    # Generate connectivity of each element based on their node numbers.
    NC = zeros(Int64, N_Elements, 4)
    for j = 1:Int(Nely / 2)
        for i = 1:Nelx
            NC[i+Nelx*(j-1), :] = [NodeNumbers[end-j+1, i], NodeNumbers[end-j+1, i+1], NodeNumbers[end-j, i+1], NodeNumbers[end-j, i]]
        end
    end

    # Degree of Freedoms of the Elements
    DOFs = reshape(1 * NC, (N_Elements, 4, 1))

    # Number of Degrees of Freedom for Helmholtz Filter
    NDOF_Filter = 1 * N_Nodes

    # Dependency Matrix
    D1 = sparse(I, N_Slaves, N_Masters)

    Ns = [N_Masters, N_Slaves, N_Internals, NDOF_Filter]

    return NodeNumbers, Ns, DOFs, D1
end

function Filter_Shape_Functions(coords, ξ, η)
    Nbar, B = zeros(Float64, 4, 2), zeros(Float64, 2, 4)

    N = 1 / 4 * [(1 - ξ) * (1 - η) (1 + ξ) * (1 - η) (1 + ξ) * (1 + η) (1 - ξ) * (1 + η)]
    N_ξ = 1 / 4 * [-(1 - η) (1 - η) (1 + η) -(1 + η)]
    N_η = 1 / 4 * [-(1 - ξ) -(1 + ξ) (1 + ξ) (1 - ξ)]

    # Construct [Nbar] (4x2)
    Nbar = [N_ξ' N_η']

    # Calculate the Jacobian
    J = (coords' * Nbar)'
    detJ, Jinv = det(J), inv(J)
    dNdX = Nbar * Jinv

    # Construct [B](2x4) Matrix
    for i = 1:2
        for j = 1:4
            B[i, j] = dNdX[j, i]
        end
    end

    return N, B, detJ

end

function Helmholtz_Filter(𝜌_Unfiltered, r)
    global F
    # ASSEMBLY
    K, F = spzeros(Float64, Ns_Filter[4], Ns_Filter[4]), zeros(Float64, Ns_Filter[4], 1)
    for el = 1:N_Elements
        Ce = vec(DOFs_Filter[el, :, :]')

        elK, elF = spzeros(Float64, 4, 4), zeros(Float64, 4, 1)
        for gp in eachrow(gps)
            ξ, η = gp[2:end]

            N, B, detJ = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)

            elK += gp[1] * (r^2 * B' * B + N' * N) * detJ
            elF += gp[1] * (𝜌_Unfiltered[el] * N') * detJ
        end
        K[Ce, Ce] += elK
        F[Ce] += elF
    end

    # SOLUTION
    Kmm = K[1:1*Ns_Filter[1], 1:1*Ns_Filter[1]]
    Kms = K[1:1*Ns_Filter[1], 1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2])]
    Kmi = K[1:1*Ns_Filter[1], 1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3])]

    Ksm = K[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2]), 1:1*Ns_Filter[1]]
    Kss = K[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2]), 1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2])]
    Ksi = K[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2]), 1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3])]

    Kim = K[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3]), 1:1*Ns_Filter[1]]
    Kis = K[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3]), 1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2])]
    Kii = K[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3]), 1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3])]

    A_Upper = hcat(Kmm + Kms * D1_Filter + D1_Filter' * Ksm + D1_Filter' * Kss * D1_Filter, Kmi + D1_Filter' * Ksi)
    A_Lower = hcat(Kim + Kis * D1_Filter, Kii)
    A_Reduced = vcat(A_Upper, A_Lower)

    Fm = F[1:1*Ns_Filter[1]]
    Fs = F[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2])]
    Fi = F[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3])]

    Bm = (Fm + D1_Filter'Fs)
    Bi = Fi
    B_Reduced = vcat(Bm, Bi)

    𝜌_Reduced = A_Reduced \ B_Reduced

    𝜌_Masters = 𝜌_Reduced[1:1*Ns_Filter[1]]
    𝜌_Slaves = D1_Filter * 𝜌_Masters
    𝜌_Internals = 𝜌_Reduced[1*Ns_Filter[1]+1:end]

    𝜌_Filtered = vcat(𝜌_Masters, 𝜌_Slaves, 𝜌_Internals)

    return 𝜌_Filtered, K
end

function constitutiveEqs(U, coords, ξ, η, 𝜌)
    λ, μ = E * v / ((1 + v) * (1 - 2 * v)), E / (2 * (1 + v))
    Nbar, B = zeros(Float64, 4, 2), zeros(Float64, 4, 8)

    N_ξ = 1 / 4 * [-(1 - η) (1 - η) (1 + η) -(1 + η)]
    N_η = 1 / 4 * [-(1 - ξ) -(1 + ξ) (1 + ξ) (1 - ξ)]

    # Construct [Nbar] (4x2)
    Nbar = [N_ξ' N_η']

    # Calculate the Jacobian
    J = coords' * Nbar
    detJ, Jinv = det(J), inv(J)
    dNdX = Nbar * Jinv

    if detJ < 0
        println("NEGATIVE DETERMINANT (detJ)")
    end

    # Construct [B](4x8) Matrix
    for i = 1:4
        for j = 1:4
            B[i, 2*(j-1)+mapI[i]] = dNdX[j, mapJ[i]]
        end
    end

    # Calculate gradients
    ∇u = reshape(U, 2, 4) * dNdX  # Displacement gradient
    F = ∇u + I                  # Deformation Gradient
    detF, Finv = det(F), inv(F)

    # First Piola-Kirchoff Stress Tensor
    P = μ * F - (μ - λ * log(detF)) * Finv'

    # Nominal stiffness matrix
    A(i, j, k, l) = μ * I[j, l] * I[i, k] + λ * Finv'[i, j] * Finv'[k, l] + (μ - λ * log(detF)) * Finv'[i, l] * Finv'[k, j]

    # Linear Stiffness Matrix
    C_Linear(i, j, k, l) = λ * I[i, j] * I[k, l] + μ * (I[i, k] * I[j, l] + I[i, l] * I[j, k])

    # Linear Stress Matrix
    S_Linear = zeros(2, 2)
    @tullio S_Linear[i, j] = C_Linear(i, j, k, l) * ∇u[k, l]

    G(i, j, k, l, m, n) = Finv'[i, j] * Finv'[k, l] * Finv'[m, n]
    # Stiffness Tensor Derivative Z
    Z(i, j, k, l, m, n) = -λ * (G(m, j, i, n, k, l) + G(i, j, m, l, k, n) + G(m, n, i, l, k, j)) - (μ - λ * log(detF)) * (G(m, l, i, n, k, j) + G(i, l, m, j, k, n))

    P_rho = g(𝜌) * (H_Steep(𝜌^p) * P + (1 - H_Steep(𝜌^p)) * S_Linear)
    A_rho(i, j, k, l) = g(𝜌) * (H_Steep(𝜌^p) * A(i, j, k, l) + (1 - H_Steep(𝜌^p)) * C_Linear(i, j, k, l))
    Z_rho(i, j, k, l, m, n) = g(𝜌) * H_Steep(𝜌^p) * Z(i, j, k, l, m, n)

    # Rewrite Tensors in Voigt Notation
    @tullio P_1D[i] := P[mapI[i], mapJ[i]]
    @tullio P_rho_1D[i] := P_rho[mapI[i], mapJ[i]]
    @tullio S_Linear_1D[i] := S_Linear[mapI[i], mapJ[i]]
    @tullio C_Linear_2D[i, j] := C_Linear(mapI[i], mapJ[i], mapI[j], mapJ[j])
    @tullio A_2D[i, j] := A(mapI[i], mapJ[i], mapI[j], mapJ[j])
    @tullio A_rho_2D[i, j] := A_rho(mapI[i], mapJ[i], mapI[j], mapJ[j])

    return B, detJ, P_1D, A_2D, P_rho_1D, A_rho_2D, S_Linear_1D, C_Linear_2D, Z_rho
end

function Displacement_Solver(𝜌_Filtered, TOL::Float64, DT::Float64)
    U_Nonlin = zeros(Float64, Ns[5], 1)
    L1, L2 = [Lx 0.0]', [0.0 Ly]'
    function Assembly(U_Nonlin::Matrix{Float64})
        # Update the elements and assemble the system matrices.
        F, K = zeros(Float64, Ns[5], 1), spzeros(Float64, Ns[5], Ns[5])
        for el = 1:N_Elements
            Ce = vec(DOFs[el, :, :]')
            Ce_Filter = vec(DOFs_Filter[el, :, :]')
            elFint, elK = zeros(Float64, 8, 1), spzeros(Float64, 8, 8)
            for gp in eachrow(gps)
                ξ, η = gp[2:end]

                N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]

                B, detJ, P_1D, A_2D, P_rho_1D, A_rho_2D, S_Linear_1D, C_Linear_2D, Z_rho = constitutiveEqs(U_Nonlin[Ce, :], elCoords[el, :, :], ξ, η, 𝜌_x)

                # Calculate the internal force matrix
                elFint += gp[1] * (B' * P_rho_1D) * detJ
                # Calculate the element tangent stiffness matrix
                elK += gp[1] * (B' * A_rho_2D * B) * detJ
            end
            F[Ce] += elFint
            K[Ce, Ce] += elK
        end
        return K, F
    end

    function solve(K::SparseMatrixCSC, F::Matrix{Float64})
        Kpp = K[1:2*Ns[1], 1:2*Ns[1]]
        Kpm = K[1:2*Ns[1], 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Kps = K[1:2*Ns[1], 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Kpi = K[1:2*Ns[1], 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

        Kmp = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 1:2*Ns[1]]
        Kmm = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Kms = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Kmi = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

        Ksp = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 1:2*Ns[1]]
        Ksm = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Kss = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Ksi = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

        Kip = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 1:2*Ns[1]]
        Kim = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Kis = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Kii = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

        K_Upper = hcat(Kpm + Kps * D1, Kpi)
        K_Middle = hcat(Kmm + Kms * D1 + D1' * Ksm + D1' * Kss * D1, Kmi + D1' * Ksi)
        K_Lower = hcat(Kim + Kis * D1, Kii)
        K_Reduced = vcat(K_Middle, K_Lower)

        Fp = F[1:2*Ns[1]]
        Fm = F[2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Fs = F[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Fi = F[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]
        R_Reduced = vcat(Fm + D1' * Fs, Fi)
        R_Norm = norm(R_Reduced)

        dU = K_Reduced \ -R_Reduced

        return dU, R_Norm
    end

    t, tFinal = (0.0, 1.0)
    A(t) = (t / tFinal) * (F_M - I)

    Rep_MU = abs(maximum(F_M - I)) * Lx

    while round(t, digits=2) < tFinal
        t += DT
        H = D2 * vcat(A(t) * L1, A(t) * L2)

        U_Prescribed = vcat(A(t) * [-Lx / 2 Ly / 2]', A(t) * [-Lx / 2 0]', A(t) * [0 Ly / 2]', A(t) * [0 0]', A(t) * [Lx / 2 Ly / 2]', A(t) * [Lx / 2 0]')
        U_Masters = U_Nonlin[2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        U_Slaves = D1 * U_Masters + H
        U_Internals = U_Nonlin[2*(Ns[1]+Ns[2]+Ns[3])+1:end]
        U_Nonlin = vcat(U_Prescribed, U_Masters, U_Slaves, U_Internals)

        K, F = Assembly(U_Nonlin)
        R_Norm = Inf
        while R_Norm / Rep_MU > TOL
            dU, R_Norm = solve(K, F)
            U_Masters += dU[1:2*Ns[2]]
            U_Internals += dU[2*Ns[2]+1:end]
            U_Slaves = D1 * U_Masters + H
            U_Nonlin = vcat(U_Prescribed, U_Masters, U_Slaves, U_Internals)
            K, F = Assembly(U_Nonlin)
        end
    end

    return U_Nonlin
end

function χ_Solver(k, l, U_Nonlin, 𝜌_Filtered)
    function Assembly()
        # Update the elements and assemble the system matrices.
        K = spzeros(Float64, Ns[5], Ns[5])
        for el = 1:N_Elements
            elK = spzeros(Float64, 8, 8)
            Ce = vec(DOFs[el, :, :]')
            Ce_Filter = vec(DOFs_Filter[el, :, :]')
            for gp in eachrow(gps)
                ξ, η = gp[2:end]

                N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]

                B, detJ, P_1D, A_2D, P_rho_1D, A_rho_2D, S_Linear_1D, C_Linear_2D, Z_rho = constitutiveEqs(U_Nonlin[Ce, :], elCoords[el, :, :], ξ, η, 𝜌_x)

                # Calculate the element tangent stiffness matrix
                elK += gp[1] * (B' * A_rho_2D * B) * detJ
            end
            K[Ce, Ce] += elK
        end
        return K
    end

    function solve(K::SparseMatrixCSC, H)
        Kpp = K[1:2*Ns[1], 1:2*Ns[1]]
        Kpm = K[1:2*Ns[1], 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Kps = K[1:2*Ns[1], 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Kpi = K[1:2*Ns[1], 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

        Kmp = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 1:2*Ns[1]]
        Kmm = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Kms = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Kmi = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

        Ksp = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 1:2*Ns[1]]
        Ksm = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Kss = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Ksi = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

        Kip = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 1:2*Ns[1]]
        Kim = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
        Kis = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
        Kii = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

        K_Upper = hcat(Kpm + Kps * D1, Kpi)
        K_Middle = hcat(Kmm + Kms * D1 + D1' * Ksm + D1' * Kss * D1, Kmi + D1' * Ksi)
        K_Lower = hcat(Kim + Kis * D1, Kii)
        K_Reduced = vcat(K_Middle, K_Lower)

        χ_Prescribed = vcat(unitStrain(k, l) * [-Lx / 2 Ly / 2]', unitStrain(k, l) * [-Lx / 2 0]', unitStrain(k, l) * [0 Ly / 2]',
            unitStrain(k, l) * [0 0]', unitStrain(k, l) * [Lx / 2 Ly / 2]', unitStrain(k, l) * [Lx / 2 0]')
        Bp = -Kps * H
        Bm = -(Kms + D1' * Kss) * H - (Kmp + D1' * Ksp) * χ_Prescribed
        Bi = -Kis * H - Kip * χ_Prescribed
        B_Reduced = vcat(Bm, Bi)

        χ_Reduced = K_Reduced \ B_Reduced

        χ_Masters = χ_Reduced[1:2*Ns[2]]
        χ_Slaves = D1 * χ_Masters + H
        χ_Internals = χ_Reduced[2*Ns[2]+1:end]

        χ = vcat(χ_Prescribed, χ_Masters, χ_Slaves, χ_Internals)
        return χ
    end

    K = Assembly()
    H = D2 * vcat(unitStrain(k, l) * e[1, :] * Lx, unitStrain(k, l) * e[2, :] * Ly)
    χ = solve(K, H)

    return χ
end

function Homogenizer_VolumeIntegral(U_Nonlin, χ, 𝜌_Filtered)
    function Homogenize_Stiffness()
        Ael = zeros(Float64, 4, 4)
        for ij in 1:4
            for kl in 1:4
                Aijkl = 0.0
                for el = 1:N_Elements
                    Ce = vec(DOFs[el, :, :]')
                    Ce_Filter = vec(DOFs_Filter[el, :, :]')
                    for gp in eachrow(gps)
                        ξ, η = gp[2:end]

                        N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                        𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]

                        B, detJ, P_1D, A_2D, P_rho_1D, A_rho_2D, S_Linear_1D, C_Linear_2D, Z_rho = constitutiveEqs(U_Nonlin[Ce, :], elCoords[el, :, :], ξ, η, 𝜌_x)

                        ∇χ⁽ⁱʲ⁾ = B * χ[ij][Ce, :]
                        ∇χ⁽ᵏˡ⁾ = B * χ[kl][Ce, :]

                        Aijkl += 2 * gp[1] * (∇χ⁽ⁱʲ⁾'*A_rho_2D*∇χ⁽ᵏˡ⁾)[1] * detJ
                    end
                end
                Ael[ij, kl] = Aijkl
            end
        end

        return 1 / (Lx * Ly) * Ael
    end

    function Homogenize_Stress()
        Pel = zeros(Float64, 4, 1)
        for el = 1:N_Elements
            Ce = vec(DOFs[el, :, :]')
            Ce_Filter = vec(DOFs_Filter[el, :, :]')
            for gp in eachrow(gps)
                ξ, η = gp[2:end]

                N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]

                B, detJ, P_1D, A_2D, P_rho_1D, A_rho_2D, S_Linear_1D, C_Linear_2D, Z_rho = constitutiveEqs(U_Nonlin[Ce, :], elCoords[el, :, :], ξ, η, 𝜌_x)

                Pel += 2 * gp[1] * (P_rho_1D) * detJ
            end

        end

        return 1 / (Lx * Ly) * Pel
    end

    return Homogenize_Stress(), Homogenize_Stiffness()
end

function Stiffness_Sensitivity_Filtered(i::Int64, j::Int64, r::Int64, s::Int64, U_Nonlin::Matrix{Float64}, χ::Vector{Matrix{Float64}}, 𝜌_Filtered)
    L1, L2 = [Lx 0.0]', [0.0 Ly]'
    function FindTrialU(ij, rs)
        function Assembly()
            K, F = spzeros(Float64, Ns[5], Ns[5]), zeros(Float64, Ns[5], 1)
            for el = 1:N_Elements
                elK, elF = spzeros(Float64, 8, 8), zeros(Float64, 8, 1)
                Ce = vec(DOFs[el, :, :]')
                Ce_Filter = vec(DOFs_Filter[el, :, :]')
                for gp in eachrow(gps)
                    ξ, η = gp[2:end]

                    N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                    𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]

                    B, detJ, P_1D, A_2D, P_rho_1D, A_rho_2D, S_Linear_1D, C_Linear_2D, Z_rho = constitutiveEqs(U_Nonlin[Ce, :], elCoords[el, :, :], ξ, η, 𝜌_x)

                    ∇χ⁽ⁱʲ⁾_Voigt = B * χ[ij][Ce, :] # In Voigt notation
                    ∇χ⁽ʳˢ⁾_Voigt = B * χ[rs][Ce, :] # In Voigt notation

                    @tullio ∇χ⁽ⁱʲ⁾[mapI[n], mapJ[n]] := ∇χ⁽ⁱʲ⁾_Voigt[n]
                    @tullio ∇χ⁽ʳˢ⁾[mapI[n], mapJ[n]] := ∇χ⁽ʳˢ⁾_Voigt[n]

                    Z_Cont = zeros(2, 2)
                    @tullio Z_Cont[ii, jj] = Z_rho(ii, jj, kk, ll, mm, nn) * ∇χ⁽ⁱʲ⁾[kk, ll] * ∇χ⁽ʳˢ⁾[mm, nn]
                    @tullio Q[n] := Z_Cont[mapI[n], mapJ[n]]

                    elF += gp[1] * (B' * Q) * detJ
                    elK += gp[1] * (B' * A_rho_2D * B) * detJ
                end
                F[Ce] += elF
                K[Ce, Ce] += elK
            end
            return K, F
        end

        function solve(K::SparseMatrixCSC, F::Matrix{Float64})
            Kpp = K[1:2*Ns[1], 1:2*Ns[1]]
            Kpm = K[1:2*Ns[1], 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
            Kps = K[1:2*Ns[1], 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
            Kpi = K[1:2*Ns[1], 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

            Kmp = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 1:2*Ns[1]]
            Kmm = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
            Kms = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
            Kmi = K[2*Ns[1]+1:2*(Ns[1]+Ns[2]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

            Ksp = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 1:2*Ns[1]]
            Ksm = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
            Kss = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
            Ksi = K[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

            Kip = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 1:2*Ns[1]]
            Kim = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*Ns[1]+1:2*(Ns[1]+Ns[2])]
            Kis = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
            Kii = K[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4]), 2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

            K_Upper = hcat(Kpm + Kps * D1, Kpi)
            K_Middle = hcat(Kmm + Kms * D1 + D1' * Ksm + D1' * Kss * D1, Kmi + D1' * Ksi)
            K_Lower = hcat(Kim + Kis * D1, Kii)
            K_Reduced = vcat(K_Middle, K_Lower)

            Fp = F[1:2*Ns[1]]
            Fm = F[2*Ns[1]+1:2*(Ns[1]+Ns[2])]
            Fs = F[2*(Ns[1]+Ns[2])+1:2*(Ns[1]+Ns[2]+Ns[3])]
            Fi = F[2*(Ns[1]+Ns[2]+Ns[3])+1:2*(Ns[1]+Ns[2]+Ns[3]+Ns[4])]

            Bm = -(Fm + D1' * Fs)
            Bi = -Fi
            B_Reduced = vcat(Bm, Bi)

            U_Reduced = K_Reduced \ B_Reduced

            U_Prescribed = zeros(2 * 6, 1)
            U_Masters = U_Reduced[1:2*Ns[2]]
            U_Slaves = D1 * U_Masters
            U_Internals = U_Reduced[2*Ns[2]+1:end]

            TrialU = vcat(U_Prescribed, U_Masters, U_Slaves, U_Internals)

            return TrialU
        end

        K, F = Assembly()
        TrialU = solve(K, F)
        return TrialU
    end

    ij, rs = convertVoigt(i, j), convertVoigt(r, s)
    TrialU = FindTrialU(ij, rs)

    A_δθ_Full_Filtered = zeros(Float64, Ns[6], 1)
    for el = 1:N_Elements
        elA_δθ = zeros(4, 1)
        Ce = vec(DOFs[el, :, :]')
        Ce_Filter = vec(DOFs_Filter[el, :, :]')
        for gp in eachrow(gps)
            ξ, η = gp[2:end]

            N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
            𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]

            B, detJ, P_1D, A_2D, P_rho_1D, A_rho_2D, S_Linear_1D, C_Linear_2D, Z_rho = constitutiveEqs(U_Nonlin[Ce, :], elCoords[el, :, :], ξ, η, 𝜌_x)

            ∇TrialU = B * TrialU[Ce, :]
            ∇χ⁽ⁱʲ⁾ = B * χ[ij][Ce, :]
            ∇χ⁽ʳˢ⁾ = B * χ[rs][Ce, :]

            Coeff1 = g_Deriv(𝜌_x) * H_Steep(𝜌_x^p) + p * 𝜌_x^(p - 1) * g(𝜌_x) * H_Steep_Deriv(𝜌_x^p)
            Coeff2 = g_Deriv(𝜌_x) - Coeff1

            δP = (Coeff1 * P_1D + Coeff2 * S_Linear_1D)
            δA = (Coeff1 * A_2D + Coeff2 * C_Linear_2D)

            elA_δθ += gp[1] * N_Filter' * (∇χ⁽ⁱʲ⁾' * δA * ∇χ⁽ʳˢ⁾ + δP' * ∇TrialU) * detJ
        end
        A_δθ_Full_Filtered[Ce_Filter] += 1 / (Lx * Ly) * elA_δθ
    end

    return A_δθ_Full_Filtered
end

function Stiffness_Sensitivity_Unfiltered(Sens, U_Nonlin, χ, 𝜌_Filtered, K_Filter)
    A_δθ_Full_Filtered = zeros(Ns[6], length(Sens))
    for n = axes(A_δθ_Full_Filtered, 2)
        i, j, r, s = Sens[n]
        A_δθ_Full_Filtered[:, n] = Stiffness_Sensitivity_Filtered(i, j, r, s, U_Nonlin, χ, 𝜌_Filtered)
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

        Fm = A_δθ_Full_Filtered[1:1*Ns_Filter[1], :]
        Fs = A_δθ_Full_Filtered[1*Ns_Filter[1]+1:1*(Ns_Filter[1]+Ns_Filter[2]), :]
        Fi = A_δθ_Full_Filtered[1*(Ns_Filter[1]+Ns_Filter[2])+1:1*(Ns_Filter[1]+Ns_Filter[2]+Ns_Filter[3]), :]

        Bm = (Fm + D1_Filter'Fs)
        Bi = Fi
        B_Reduced = vcat(Bm, Bi)

        w_Reduced = A_Reduced \ B_Reduced

        w_Masters = w_Reduced[1:1*Ns_Filter[1], :]
        w_Slaves = D1_Filter * w_Masters
        w_Internals = w_Reduced[1*Ns_Filter[1]+1:end, :]

        w_Filter = vcat(w_Masters, w_Slaves, w_Internals)
    end

    A_δθ_Full_Unfiltered = zeros(Float64, N_Elements, size(A_δθ_Full_Filtered, 2))
    for i = axes(A_δθ_Full_Unfiltered, 2)
        for el = 1:N_Elements
            el_δθ = 0.0
            Ce_Filter = vec(DOFs_Filter[el, :, :]')
            for gp in eachrow(gps)
                ξ, η = gp[2:end]
                N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                el_δθ += gp[1] * (N_Filter*w_Filter[:, i][Ce_Filter])[1] * detJ_Filter
            end
            A_δθ_Full_Unfiltered[el, i] = el_δθ
        end
    end

    A_δθ_Unfiltered = zeros(Int(N_Elements / 2), size(A_δθ_Full_Unfiltered, 2))
    elements_BL = reduce(hcat, [vcat(Int(i * Nelx + 1):Int(i * Nelx + Nelx / 2)) for i = 0:Nely/2-1])'
    for l = axes(A_δθ_Unfiltered, 2)
        for k = 1:Int(N_Elements / 2)
            el = elements_BL[k]
            i, j = Int(ceil(el / Nely)), (el - 1) % Nely + 1
            elV = el + Nelx - (2 * j - 1)
            A_δθ_Unfiltered[k, l] = A_δθ_Full_Unfiltered[el, l] + A_δθ_Full_Unfiltered[elV, l]
        end
    end
    return eachcol(2 * A_δθ_Unfiltered)
end

#= function Stress_Sensitivity_Filtered(i, j, U_Nonlin, χ, 𝜌_Filtered)
    ij = convertVoigt(i, j)
    P_δθ_Full_Filtered = zeros(Float64, Ns[5], 1)
    for el = 1:N_Elements
        el_δθ = zeros(4, 1)
        Ce = vec(DOFs[el, :, :]')
        Ce_Filter = vec(DOFs_Filter[el, :, :]')
        for gp in eachrow(gps)
            ξ, η = gp[2:end]

            N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
            𝜌_x = (N_Filter*𝜌_Filtered[Ce_Filter])[1]

            B, detJ, P_1D, A_2D, P_rho_1D, A_rho_2D, S_Linear_1D, C_Linear_2D, Z_rho = constitutiveEqs(U_Nonlin[Ce, :], elCoords[el, :, :], ξ, η, 𝜌_x)

            ∇χ⁽ⁱʲ⁾ = B * χ[ij][Ce, :]

            Coeff1 = g_Deriv(𝜌_x) * H_Steep(𝜌_x^p) + p * 𝜌_x^(p - 1) * g(𝜌_x) * H_Steep_Deriv(𝜌_x^p)
            Coeff2 = g_Deriv(𝜌_x) - Coeff1

            δP = (Coeff1 * P_1D + Coeff2 * S_Linear_1D)

            el_δθ += gp[1] * N_Filter' * (δP' * ∇χ⁽ⁱʲ⁾) * detJ
        end
        P_δθ_Full_Filtered[Ce_Filter] += 1 / (Lx * Ly) * el_δθ
    end

    return P_δθ_Full_Filtered
end

function Stress_Sensitivity_Unfiltered(Sens, U_Nonlin, χ, 𝜌_Filtered, K_Filter)
    P_δθ_Full_Filtered = zeros(Ns[5], length(Sens))
    for n = axes(P_δθ_Full_Filtered, 2)
        i, j = Sens[n]
        P_δθ_Full_Filtered[:, n] = Stress_Sensitivity_Filtered(i, j, U_Nonlin, χ, 𝜌_Filtered)
    end

    w_Filter = K_Filter \ P_δθ_Full_Filtered

    P_δθ_Full_Unfiltered = zeros(Float64, N_Elements, size(P_δθ_Full_Filtered, 2))
    for i = axes(P_δθ_Full_Unfiltered, 2)
        for el = 1:N_Elements
            elP_δθ = 0.0
            Ce_Filter = vec(DOFs_Filter[el, :, :]')
            for gp in eachrow(gps)
                ξ, η = gp[2:end]
                N_Filter, B_Filter, detJ_Filter = Filter_Shape_Functions(elCoords[el, :, :], ξ, η)
                elP_δθ += gp[1] * (N_Filter*w_Filter[:, i][Ce_Filter])[1] * detJ_Filter
            end
            P_δθ_Full_Unfiltered[el, i] = elP_δθ
        end
    end

    elements_BL = reduce(hcat, [vcat(Int(i * Nelx + 1):Int(i * Nelx + Nelx / 2)) for i = 0:Nely/2-1])'
    P_δθ_Unfiltered = zeros(Int(N_Elements / 4), size(P_δθ_Full_Unfiltered, 2))
    for l = axes(P_δθ_Unfiltered, 2)
        for k = 1:Int(N_Elements / 4)
            el = elements_BL[k]
            i, j = Int(ceil(el / Nely)), (el - 1) % Nely + 1
            elV, elH, elVH = i * Nelx - j + 1, (Nely - i) * Nelx + j, (Nely - i) * Nelx + Nelx - j + 1
            P_δθ_Unfiltered[k, l] = P_δθ_Full_Unfiltered[el, l] + P_δθ_Full_Unfiltered[elV, l] + P_δθ_Full_Unfiltered[elH, l] + P_δθ_Full_Unfiltered[elVH, l]
        end
    end
    return eachcol(P_δθ_Unfiltered)
end =#