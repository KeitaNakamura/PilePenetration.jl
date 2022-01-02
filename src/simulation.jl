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
    friction_with_pile::Float64
end

struct PointState
    m::Float64
    V::Float64
    x::Vec{2, Float64}
    v::Vec{2, Float64}
    b::Vec{2, Float64}
    σ::SymmetricSecondOrderTensor{3, Float64, 6}
    ϵ::SymmetricSecondOrderTensor{3, Float64, 6}
    ∇v::SecondOrderTensor{3, Float64, 9}
    C::Mat{2, 3, Float64, 6}
    side_length::Vec{2, Float64}
    friction_with_pile_inside::Float64
    friction_with_pile_outside::Float64
    index::Int
    layerindex::Int
end


function main_simulation()::Cint
    inputtoml = only(ARGS)
    try
        main_simulation(inputtoml)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function main_simulation(inputtoml::AbstractString)
    proj_dir = splitdir(inputtoml)[1]
    INPUT = parseinput(inputtoml)
    output_dir = joinpath(proj_dir, INPUT.General.output_folder_name)
    mkpath(output_dir)
    cp(inputtoml, joinpath(output_dir, "input.toml"); force = true)
    main_simulation(proj_dir, INPUT)
end

function main_simulation(proj_dir::AbstractString, INPUT::NamedTuple)

    # General
    (xmin, xmax), (ymin, ymax) = INPUT.General.domain
    dx = INPUT.General.grid_space
    g = INPUT.General.gravity
    total_time = INPUT.General.total_time

    # SoilLayer
    soillayers = INPUT.SoilLayer # reorder layers from low to high
    H = sum(layer -> layer.thickness, soillayers)
    @assert H ≤ ymax

    # Pile
    D_i = INPUT.Pile.diameter_head
    d_i = INPUT.Pile.diameter_tip
    pile_length = INPUT.Pile.pile_length
    tapered_length = INPUT.Pile.tapered_length
    thickness = INPUT.Pile.thickness
    vy_pile = INPUT.Pile.velocity
    vacuum = INPUT.Pile.vacuum
    vacuum_height = INPUT.Pile.vacuum_height
    @assert thickness ≥ dx

    # Advanced
    CFL = INPUT.Advanced.CFL
    contact_threshold_scale = INPUT.Advanced.contact_threshold_scale
    contact_penalty_parameter = INPUT.Advanced.contact_penalty_parameter


    grid = Grid(NodeState, LinearWLS(QuadraticBSpline()), xmin:dx:xmax, ymin:dx:ymax; coordinate_system = Axisymmetric())
    pointstate = generate_pointstate((x,y) -> y < H, PointState, grid)
    cache = MPCache(grid, pointstate.x)

    layermodels = map(soillayers) do layer
        E = layer.youngs_modulus
        ν = layer.poissons_ratio
        ϕ = layer.friction_angle
        ψ = layer.dilatancy_angle
        c = layer.cohesion
        tension_cutoff = layer.tension_cutoff
        elastic = LinearElastic(; E, ν)
        DruckerPrager(elastic, :circumscribed; c, ϕ, ψ, tension_cutoff)
    end
    bottom = ymin
    for i in length(soillayers):-1:1 # from low to high
        layer = soillayers[i]
        Threads.@threads for p in 1:length(pointstate)
            y = pointstate.x[p][2]
            if bottom ≤ y ≤ bottom + layer.thickness
                pointstate.layerindex[p] = i
            end
        end
        bottom += layer.thickness
    end

    Threads.@threads for p in 1:length(pointstate)
        layerindex = pointstate.layerindex[p]
        σ_y = 0.0
        for layer in soillayers[begin:layerindex-1]
            γ₀ = layer.unit_weight
            σ_y += -γ₀ * layer.thickness
        end
        h = sum(layer -> layer.thickness, soillayers[layerindex:end])
        layer = soillayers[layerindex]
        y = pointstate.x[p][2]
        γ0 = layer.unit_weight
        ν = layer.poissons_ratio
        σ_y += -γ0 * (h - y)
        σ_x = σ_y * ν / (1 - ν)
        pointstate.σ[p] = (@Mat [σ_x 0.0 0.0
                                 0.0 σ_y 0.0
                                 0.0 0.0 σ_x]) |> symmetric
        pointstate.m[p] = (γ0/g) * pointstate.V[p]
        pointstate.friction_with_pile_inside[p] = layer.friction_with_pile_inside
        pointstate.friction_with_pile_outside[p] = layer.friction_with_pile_outside
    end
    @. pointstate.b = Vec(0.0, -g)
    Poingr.reorder_pointstate!(pointstate, cache)

    R_i = D_i / 2 # radius
    r_i = d_i / 2 # radius
    R_o = R_i + thickness
    r_o = r_i + thickness
    pile = Polygon([Vec(R_i, pile_length), Vec(R_i, tapered_length), Vec(r_i, 0.0),
                    Vec(r_o, 0.0), Vec(R_o, tapered_length), Vec(R_o, pile_length)])
    pile_inside = Polygon([Vec(0.0, pile_length), Vec(0.0, 0.0),
                           Vec(r_i, 0.0), Vec(R_i, tapered_length), Vec(R_i, pile_length)])
    translate!(pile, Vec(0.0, H+dx))
    translate!(pile_inside, Vec(0.0, H+dx))
    v_pile = Vec(0.0, -vy_pile)

    pile_center_0 = centroid(pile)
    ground_height_0 = find_ground_pos(pointstate.x, gridsteps(grid, 1))

    # Output files
    outputs = Dict{String, Any}()
    output_dir = joinpath(proj_dir, INPUT.General.output_folder_name)
    outputs["output directory"] = output_dir
    ## serialize
    mkpath(joinpath(output_dir, "serialize"))
    ## paraview
    mkpath(joinpath(output_dir, "paraview"))
    outputs["paraview file"] = joinpath(output_dir, "paraview", "output")
    paraview_collection(vtk_save, outputs["paraview file"])
    ## history
    outputs["history file"] = joinpath(output_dir, "history.csv")
    outputhistory_head(outputs["history file"])

    println("Start: ", now())
    println("Particles: ", length(pointstate))

    t = 0.0
    logger = Logger(0.0:INPUT.General.output_interval:total_time; progress = INPUT.General.show_progress)
    update!(logger, t)
    writeoutput(outputs, grid, pointstate, pile, logger, ground_height_0, pile_center_0, t, INPUT)
    while !isfinised(logger, t)
        dt = minimum(pointstate) do p
            ρ = p.m / p.V
            elastic = layermodels[p.layerindex].elastic
            vc = soundspeed(elastic.K, elastic.G, ρ)
            CFL * dx / vc
        end

        update!(cache, grid, pointstate.x)
        P2G!(grid, pointstate, cache, pile, v_pile, dt, INPUT)
        for bd in eachboundary(grid)
            @inbounds grid.state.v[bd.I] = boundary_velocity(grid.state.v[bd.I], bd.n)
        end
        G2P!(pointstate, grid, cache, layermodels, pile, dt)

        translate!(pile, v_pile * dt)
        translate!(pile_inside, v_pile * dt)
        update!(logger, t += dt)

        if islogpoint(logger)
            Poingr.reorder_pointstate!(pointstate, cache)
            if vacuum
                inds = findall(xₚ -> (xₚ in pile_inside) && xₚ[2] > (tip_height(pile) + vacuum_height), pointstate.x)
                deleteat!(pointstate, inds)
            end
            writeoutput(outputs, grid, pointstate, pile, logger, ground_height_0, pile_center_0, t, INPUT)
        end
    end
end

tip_height(pile) = pile[3][2]
find_ground_pos(xₚ, dx) = maximum(x -> x[2], filter(x -> x[1] < dx, xₚ))

function P2G!(grid::Grid, pointstate::AbstractVector, cache::MPCache, pile::Polygon, v_pile, dt, INPUT)
    contact_threshold_scale = INPUT.Advanced.contact_threshold_scale
    contact_penalty_parameter = INPUT.Advanced.contact_penalty_parameter

    default_point_to_grid!(grid, pointstate, cache)

    mask = @. distance($Ref(pile), pointstate.x, minimum(pointstate.side_length)/2 * contact_threshold_scale) !== nothing
    point_to_grid!((grid.state.fcn, grid.state.vᵣ, grid.state.w_pile, grid.state.friction_with_pile), cache, mask) do it, p, i
        @_inline_meta
        @_propagate_inbounds_meta
        N = it.N
        w = it.w
        vₚ = pointstate.v[p]
        fcn = N * contact_normal_force(pile, pointstate.x[p], pointstate.m[p], pointstate.side_length[p] * contact_threshold_scale, contact_penalty_parameter, dt)
        vᵣ = w * (vₚ - v_pile)
        if fcn[1] > 0
            μ = w * pointstate.friction_with_pile_inside[p]
        else
            μ = w * pointstate.friction_with_pile_outside[p]
        end
        fcn, vᵣ, w, μ
    end
    @dot_threads grid.state.vᵣ /= grid.state.w_pile
    @dot_threads grid.state.friction_with_pile /= grid.state.w_pile
    @dot_threads grid.state.fc = contact_force(grid.state.vᵣ, grid.state.fcn, grid.state.m, dt, grid.state.friction_with_pile)

    @dot_threads grid.state.v += ((grid.state.f + grid.state.fc) / grid.state.m) * dt
end

function G2P!(pointstate::AbstractVector, grid::Grid, cache::MPCache, layermodels::Vector{<: DruckerPrager}, pile::Polygon, dt)
    default_grid_to_point!(pointstate, grid, cache, dt)
    @inbounds Threads.@threads for p in eachindex(pointstate)
        model = layermodels[pointstate.layerindex[p]]
        ∇v = pointstate.∇v[p]
        σ_n = pointstate.σ[p]
        dϵ = symmetric(∇v) * dt
        σ = update_stress(model, σ_n, dϵ)
        σ = Poingr.jaumann_stress(σ, σ_n, ∇v, dt)
        if mean(σ) > model.tension_cutoff
            # In this case, since the soil particles are not contacted with
            # each other, soils should not act as continuum.
            # This means that the deformation based on the contitutitive model
            # no longer occurs.
            # So, in this process, we just calculate the elastic strain to keep
            # the consistency with the stress which is on the edge of the yield
            # function, and ignore the plastic strain to prevent excessive generation.
            # If we include this plastic strain, the volume of the material points
            # will continue to increase unexpectedly.
            σ_tr = update_stress(model.elastic, σ_n, dϵ)
            σ = Poingr.tension_cutoff(model, σ_tr)
            dϵ = model.elastic.Dinv ⊡ (σ - σ_n)
        end
        pointstate.σ[p] = σ
        pointstate.ϵ[p] += dϵ
        pointstate.V[p] *= exp(tr(dϵ))
    end
end

function writeoutput(outputs::Dict{String, Any}, grid::Grid, pointstate::AbstractVector, pile::Polygon, logger, ground_height_0::Real, pile_center_0::Vec, t::Real, INPUT)
    output_dir = outputs["output directory"]
    paraview_file = outputs["paraview file"]

    paraview_collection(paraview_file, append = true) do pvd
        vtk_multiblock(string(paraview_file, logindex(logger))) do vtm
            vtk_points(vtm, pointstate.x) do vtk
                ϵ = pointstate.ϵ
                vtk["velocity"] = pointstate.v
                vtk["mean stress"] = @dot_lazy -mean(pointstate.σ)
                vtk["deviatoric stress"] = @dot_lazy deviatoric_stress(pointstate.σ)
                vtk["volumetric strain"] = @dot_lazy volumetric_strain(ϵ)
                vtk["deviatoric strain"] = @dot_lazy deviatoric_strain(ϵ)
                vtk["stress"] = @dot_lazy -pointstate.σ
                vtk["strain"] = ϵ
                vtk["density"] = @dot_lazy pointstate.m / pointstate.V
                vtk["soil layer"] = pointstate.layerindex
            end
            vtk_grid(vtm, pile)
            vtk_grid(vtm, grid) do vtk
                vtk["nodal force"] = vec(grid.state.f)
                vtk["nodal contact force"] = vec(grid.state.fc)
            end
            pvd[t] = vtm
        end
    end

    outputhistory_append(
        outputs["history file"],
        grid,
        pointstate,
        pile,
        INPUT.Pile.diameter_head,
        INPUT.Pile.tapered_length,
        ground_height_0,
        pile_center_0,
    )

    serialize(joinpath(output_dir, "serialize", string("save", logindex(logger))),
              (; pointstate, grid, pile, t))
end

function contact_normal_force(poly::Polygon, x::Vec{dim, T}, m::T, l::Vec{dim, T}, ξ, dt::T) where {dim, T}
    thresh = minimum(l) / 2
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
