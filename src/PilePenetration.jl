module PilePenetration

using Poingr, GeometricObjects
using StructArrays
using DelimitedFiles
using Serialization
using TOML
using Dates

using Base: @_propagate_inbounds_meta, @_inline_meta

struct NodeState
    f::Vec{2, Float64}
    fc::Vec{2, Float64}
    fcn::Vec{2, Float64}
    w::Float64
    m::Float64
    v::Vec{2, Float64}
    vᵣ::Vec{2, Float64}
    w_pile::Float64
end

struct PointState
    m::Float64
    V0::Float64
    x::Vec{2, Float64}
    v::Vec{2, Float64}
    b::Vec{2, Float64}
    σ::SymmetricSecondOrderTensor{3, Float64, 6}
    σ0::SymmetricSecondOrderTensor{3, Float64, 6}
    F::SecondOrderTensor{3, Float64, 9}
    ∇v::SecondOrderTensor{3, Float64, 9}
    C::Mat{2, 3, Float64, 6}
    side_length::Vec{2, Float64}
    index::Int
end

function julia_main()::Cint
    isempty(ARGS) && throw(ArgumentError("input.toml must be given as the first argument"))
    inputtoml = ARGS[1]
    try
        main(inputtoml)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function main(inputtoml::AbstractString)
    dict = TOML.parsefile(inputtoml)
    list = map(collect(keys(dict))) do section
        subdict = dict[section]
        Symbol(section) => (; (Symbol(key) => value for (key, value) in subdict)...)
    end
    input = (; list...)
    main(splitdir(inputtoml)[1], input)
end

function main(proj_dir::AbstractString, INPUT::NamedTuple)

    # General
    (xmin, xmax), (ymin, ymax) = INPUT.General.domain
    dx = INPUT.General.grid_space
    g = INPUT.General.gravity
    total_time = INPUT.General.total_time

    # Ground
    h = INPUT.Ground.height
    ρ₀ = INPUT.Ground.dry_density
    E = INPUT.Ground.youngs_modulus
    ν = INPUT.Ground.poissons_ratio
    ϕ = INPUT.Ground.friction_angle
    ψ = INPUT.Ground.dilatancy_angle

    # Pile
    D_i = INPUT.Pile.diameter_head
    d_i = INPUT.Pile.diameter_tip
    pile_length = INPUT.Pile.pile_length
    tapered_length = INPUT.Pile.tapered_length
    thickness = INPUT.Pile.thickness
    μ = INPUT.Pile.friction_with_ground
    vy_pile = INPUT.Pile.velocity
    vacuum = INPUT.Pile.vacuum
    vacuum_height = INPUT.Pile.vacuum_height
    @assert thickness ≥ dx

    # Advanced
    CFL = INPUT.Advanced.CFL
    contact_threshold_scale = INPUT.Advanced.contact_threshold_scale
    contact_penalty_parameter = INPUT.Advanced.contact_penalty_parameter


    coord_system = Axisymmetric()
    grid = Grid(NodeState, LinearWLS(QuadraticBSpline()), xmin:dx:xmax, ymin:dx:ymax)
    pointstate = generate_pointstate((x,y) -> y < h, PointState, grid, coord_system)
    cache = MPCache(grid, pointstate.x)
    elastic = LinearElastic(; E, ν)
    model = DruckerPrager(elastic, :circumscribed; c = 0.0, ϕ, ψ)

    for p in 1:length(pointstate)
        y = pointstate.x[p][2]
        σ_y = -ρ₀ * g * (h - y)
        σ_x = σ_y * ν / (1 - ν)
        pointstate.σ[p] = (@Mat [σ_x 0.0 0.0
                                 0.0 σ_y 0.0
                                 0.0 0.0 σ_x]) |> symmetric
    end
    @. pointstate.m = ρ₀ * pointstate.V0
    @. pointstate.F = one(SecondOrderTensor{3,Float64})
    @. pointstate.b = Vec(0.0, -g)
    @. pointstate.σ0 = pointstate.σ
    Poingr.reorder_pointstate!(pointstate, cache)

    R_i = D_i / 2 # radius
    r_i = d_i / 2 # radius
    R_o = R_i + thickness
    r_o = r_i + thickness
    tip_height_0 = h + dx
    pile = Polygon([Vec(R_i, tip_height_0 + pile_length),
                    Vec(R_i, tip_height_0 + tapered_length),
                    Vec(r_i, tip_height_0),
                    Vec(r_o, tip_height_0),
                    Vec(R_o, tip_height_0 + tapered_length),
                    Vec(R_o, tip_height_0 + pile_length)])
    inside_pile = Polygon([Vec(0.0, tip_height_0 + pile_length),
                           Vec(0.0, tip_height_0),
                           Vec(r_i, tip_height_0),
                           Vec(R_i, tip_height_0 + tapered_length),
                           Vec(R_i, tip_height_0 + pile_length)])
    v_pile = Vec(0.0, -vy_pile)

    pile_center_0 = centroid(pile)

    find_ground_pos(xₚ) = maximum(x -> x[2], filter(x -> x[1] < dx, xₚ))
    ground_pos0 = find_ground_pos(pointstate.x)

    # Output files
    ## proj
    output_dir = joinpath(proj_dir, INPUT.General.output_folder_name)
    mkpath(output_dir)

    ## paraview
    paraview_file = joinpath(output_dir, "out")
    paraview_collection(vtk_save, paraview_file)

    ## history
    csv_file = joinpath(output_dir, "history.csv")
    open(csv_file, "w") do io
        writedlm(io, ["disp" "force" "disp_inside_pile" "tip_resistance" "inside_resistance" "outside_resistance"], ',')
    end

    ## forces
    mkpath(joinpath(output_dir, "force_tip"))
    mkpath(joinpath(output_dir, "force_inside"))
    mkpath(joinpath(output_dir, "force_outside"))

    logger = Logger(0.0:INPUT.General.output_interval:total_time; progress = true)
    t = 0.0

    println("Start: ", now())
    println("Particles: ", length(pointstate))

    while !isfinised(logger, t)

        dt = minimum(pointstate) do p
            ρ = p.m / (p.V0 * det(p.F))
            vc = soundspeed(elastic.K, elastic.G, ρ)
            CFL * dx / vc
        end

        update!(cache, grid, pointstate.x)
        P2G!(grid, pointstate, cache, pile, v_pile, dt, μ, contact_threshold_scale, contact_penalty_parameter, coord_system)
        for bd in eachboundary(grid)
            @. grid.state.v[bd.indices] = boundary_velocity(grid.state.v[bd.indices], bd.n)
        end
        G2P!(pointstate, grid, cache, model, pile, dt, coord_system)

        translate!(pile, v_pile * dt)
        translate!(inside_pile, v_pile * dt)
        update!(logger, t += dt)

        if islogpoint(logger)
            Poingr.reorder_pointstate!(pointstate, cache)

            if vacuum
                tip_height = pile[3][2]
                inds = findall(xₚ -> (xₚ in inside_pile) && xₚ[2] > (tip_height + vacuum_height), pointstate.x)
                StructArrays.foreachfield(v -> deleteat!(v, inds), pointstate)
            end

            paraview_collection(paraview_file, append = true) do pvd
                vtk_multiblock(string(paraview_file, logindex(logger))) do vtm
                    vtk_points(vtm, pointstate.x) do vtk
                        ϵₚ = @dot_lazy symmetric(pointstate.F - $Ref(I))
                        vtk["velocity"] = pointstate.v
                        vtk["mean stress"] = @dot_lazy -mean(pointstate.σ)
                        vtk["deviatoric stress"] = @dot_lazy deviatoric_stress(pointstate.σ)
                        vtk["volumetric strain"] = @dot_lazy volumetric_strain(ϵₚ)
                        vtk["deviatoric strain"] = @dot_lazy deviatoric_strain(ϵₚ)
                        vtk["stress"] = @dot_lazy -pointstate.σ
                        vtk["strain"] = ϵₚ
                        vtk["density"] = @dot_lazy pointstate.m / (pointstate.V0 * det(pointstate.F))
                    end
                    vtk_grid(vtm, pile)
                    vtk_grid(vtm, grid) do vtk
                        vtk["nodal force"] = vec(grid.state.f)
                        vtk["nodal contact force"] = vec(grid.state.fc)
                    end
                    pvd[t] = vtm
                end
            end

            tip, inside, outside = extract_contact_forces(grid.state.fc, grid, pile)
            open(csv_file, "a") do io
                disp = norm(centroid(pile) - pile_center_0)
                force = -sum(grid.state.fc)[2] * 2π
                disp_inside_pile = -(find_ground_pos(pointstate.x) - ground_pos0)
                writedlm(io, [disp force disp_inside_pile sum(@view tip[:,3]) sum(@view inside[:,3]) sum(@view outside[:,3])], ',')
            end
            open(io -> writedlm(io, tip, ','), joinpath(output_dir, "force_tip", "force_tip_$(logindex(logger)).csv"), "w")
            open(io -> writedlm(io, inside, ','), joinpath(output_dir, "force_inside", "force_inside_$(logindex(logger)).csv"), "w")
            open(io -> writedlm(io, outside, ','), joinpath(output_dir, "force_outside", "force_outside_$(logindex(logger)).csv"), "w")

            serialize(joinpath(output_dir, string("save", logindex(logger))),
                      Dict("pointstate" => pointstate,
                           "grid" => grid,
                           "pile" => pile))
        end
    end
end

function P2G!(grid::Grid, pointstate::AbstractVector, cache::MPCache, pile::Polygon, v_pile, dt, μ, contact_threshold_scale, contact_penalty_parameter, coord_system)
    default_point_to_grid!(grid, pointstate, cache, coord_system)
    @dot_threads grid.state.v += (grid.state.f / grid.state.m) * dt
    mask = @. distanceto($Ref(pile), pointstate.x, pointstate.side_length * contact_threshold_scale) !== nothing
    point_to_grid!((grid.state.fcn, grid.state.vᵣ, grid.state.w_pile), cache, mask) do it, p, i
        @_inline_meta
        @_propagate_inbounds_meta
        N = it.N
        w = it.w
        vₚ = pointstate.v[p]
        fcnₚ = contact_normal_force(pile, pointstate.x[p], pointstate.m[p], pointstate.side_length[p] * contact_threshold_scale, contact_penalty_parameter, dt)
        fcn = N * fcnₚ
        vᵣ = w * (vₚ - v_pile)
        fcn, vᵣ, w
    end
    @dot_threads grid.state.vᵣ /= grid.state.w_pile
    @dot_threads grid.state.fc = contact_force(grid.state.vᵣ, grid.state.fcn, grid.state.m, dt, μ)
    @dot_threads grid.state.v += (grid.state.fc / grid.state.m) * dt
end

function G2P!(pointstate::AbstractVector, grid::Grid, cache::MPCache, model::DruckerPrager, pile::Polygon, dt, coord_system)
    Dinv = inv(model.elastic.D)
    default_grid_to_point!(pointstate, grid, cache, dt, coord_system)
    @inbounds Threads.@threads for p in eachindex(pointstate)
        σ = update_stress(model, pointstate.σ[p], symmetric(pointstate.∇v[p]) * dt)
        σ = Poingr.jaumann_stress(σ, pointstate.σ[p], pointstate.∇v[p], dt)
        F = pointstate.F[p]
        if mean(σ) > 0
            σ = zero(σ)
            ϵv = tr(Dinv ⊡ (σ - pointstate.σ0[p]))
            J = exp(ϵv)
            F = J^(1/3) * one(F)
        end
        if det(F) < 0
            F = zero(F)
        end
        pointstate.σ[p] = σ
        pointstate.F[p] = F
    end
end

function extract_contact_forces(fcᵢ, grid, pile)
    inside = Float64[]
    outside = Float64[]
    tip = Float64[]
    tip_height = pile[3][2]
    for I in eachindex(fcᵢ) # walk from lower height
        x = grid[I][1]
        y = grid[I][2] - tip_height
        fcy = -2π * fcᵢ[I][2]
        iszero(fcy) && continue
        if y < gridsteps(grid, 2)
            push!(tip, x)
            push!(tip, y)
            push!(tip, fcy)
        else
            line1 = Line(((GeometricObjects.getline(pile, 1) + reverse(GeometricObjects.getline(pile, 5))) / 2)...)
            line2 = Line(((GeometricObjects.getline(pile, 2) + reverse(GeometricObjects.getline(pile, 4))) / 2)...)
            inner = GeometricObjects.ray_casting_to_right(line1, grid[I]) ||
                GeometricObjects.ray_casting_to_right(line2, grid[I])
            if inner
                push!(inside, x)
                push!(inside, y)
                push!(inside, fcy)
            else
                push!(outside, x)
                push!(outside, y)
                push!(outside, fcy)
            end
        end
    end
    reshape_data = V -> reshape(V, 3, length(V)÷3)'
    map(reshape_data, (tip, inside, outside))
end

function distanceto(poly::Polygon, x::Vec{dim}, l::Vec{dim}) where {dim}
    thresh = mean(l) / 2
    distance(poly, x, thresh)
end

function contact_normal_force(poly::Polygon, x::Vec{dim, T}, m::T, l::Vec{dim, T}, ξ, dt::T) where {dim, T}
    thresh = mean(l) / 2
    d = distance(poly, x, thresh)
    d === nothing && return zero(Vec{dim, T})
    norm_d = norm(d)
    n = d / norm_d
    (1-ξ) * 2m/dt^2 * (thresh - norm_d) * n
end

function contact_force(vᵣ::Vec, f_n::Vec, m::Real, dt::Real, μ::Real)
    iszero(f_n) && return zero(f_n)
    n = f_n / norm(f_n)
    vᵣ_t = vᵣ - (vᵣ ⋅ n) * n
    f_t = (m / dt) * vᵣ_t
    Contact(:friction, μ, sep = true)(f_n + f_t, n)
end

function boundary_velocity(v::Vec, n::Vec)
    if n == Vec(0, -1) # bottom
        v + Contact(:sticky)(v, n)
    else
        v + Contact(:slip)(v, n)
    end
end

#=
using Plots, DelimitedFiles
plot((arr = readdlm("output.tmp/history.csv", ',', skipstart=1); @show size(arr, 1); (arr[:,1], arr[:,[2,4,5,6]]))..., label = ["total" "tip" "inside" "outside"], xlims = (0,2), ylims = (0,60e3))
=#

end # module
