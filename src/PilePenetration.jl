module PilePenetration

using Poingr
using GeometricObjects
using GeometricObjects: getline

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

function julia_main()::Cint
    if isempty(ARGS)
        inputtoml = "input.toml"
    else
        inputtoml = ARGS[1]
    end
    try
        main(inputtoml)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function parseinput(dict::Dict)
    dict2namedtuple(x::Dict) = (; (Symbol(key) => value for (key, value) in x)...)
    list = map(collect(keys(dict))) do section
        content = dict[section]
        if content isa Dict
            Symbol(section) => dict2namedtuple(content)
        elseif content isa Vector
            Symbol(section) => map(dict2namedtuple, content)
        else
            error("unreachable")
        end
    end
    (; list...)
end
function parseinput(inputtoml::AbstractString)
    parseinput(TOML.parsefile(inputtoml))
end

function main(inputtoml::AbstractString)
    proj_dir = splitdir(inputtoml)[1]
    INPUT = parseinput(inputtoml)
    output_dir = joinpath(proj_dir, INPUT.General.output_folder_name)
    mkpath(output_dir)
    cp(inputtoml, joinpath(output_dir, "input.toml"); force = true)
    main(proj_dir, INPUT)
end

function main(proj_dir::AbstractString, INPUT::NamedTuple)

    # General
    coordinate_system = INPUT.General.coordinate_system
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


    grid = Grid(NodeState, LinearWLS(QuadraticBSpline()), xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
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
    open(outputs["history file"], "w") do io
        writedlm(io, ["disp" "force" "disp_inside_pile" "tip" "inside" "outside" "taper" "straight" "tip_design" "inside_design" "outside_design" "taper_design" "straight_design"], ',')
    end
    ## forces
    outputs["force inside directory"] = joinpath(output_dir, "force_inside")
    outputs["force outside directory"] = joinpath(output_dir, "force_outside")
    mkpath(outputs["force inside directory"])
    mkpath(outputs["force outside directory"])

    println("Start: ", now())
    println("Particles: ", length(pointstate))

    t = 0.0
    logger = Logger(0.0:INPUT.General.output_interval:total_time; progress = INPUT.General.show_progress)
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
                StructArrays.foreachfield(v -> deleteat!(v, inds), pointstate)
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

    mask = @. distanceto($Ref(pile), pointstate.x, pointstate.side_length * contact_threshold_scale) !== nothing
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
    history_file = outputs["history file"]
    force_inside_directory = outputs["force inside directory"]
    force_outside_directory = outputs["force outside directory"]

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

    inside_total, outside_total = extract_contact_forces(grid.state.fc, grid, pile)

    open(history_file, "a") do io
        disp = norm(centroid(pile) - pile_center_0)
        force = -sum(grid.state.fc)[2] * 2π
        disp_inside_pile = -(find_ground_pos(pointstate.x, gridsteps(grid, 1)) - ground_height_0)
        tip, inside, outside = divide_force_into_tip_inside_outside(gridsteps(grid, 2), inside_total, outside_total)
        tip_design, inside_design, outside_design = divide_force_into_tip_inside_outside(INPUT.Pile.diameter_head, inside_total, outside_total)
        tip_taper, _, _ = divide_force_into_tip_inside_outside(INPUT.Pile.tapered_length, inside_total, outside_total)
        taper = tip_taper - tip
        straight = force - tip_taper
        taper_design = tip_taper - tip_design
        straight_design = force - tip_taper
        writedlm(io, [disp force disp_inside_pile tip inside outside taper straight tip_design inside_design outside_design taper_design straight_design], ',')
    end
    open(joinpath(force_inside_directory, "force_inside_$(logindex(logger)).csv"), "w") do io
        writedlm(io, ["height" "force"], ',')
        writedlm(io, inside_total, ',')
    end
    open(joinpath(force_outside_directory, "force_outside_$(logindex(logger)).csv"), "w") do io
        writedlm(io, ["height" "force"], ',')
        writedlm(io, outside_total, ',')
    end

    serialize(joinpath(output_dir, "serialize", string("save", logindex(logger))),
              (; pointstate, grid, pile))
end

function divide_force_into_tip_inside_outside(height_thresh, inside_total, outside_total)
    index_inside = searchsortedfirst(view(inside_total, :, 1), height_thresh) - 1
    index_outside = searchsortedfirst(view(outside_total, :, 1), height_thresh) - 1
    tip = sum(@view inside_total[begin:index_inside, 2]) + sum(@view outside_total[begin:index_outside, 2])
    inside = sum(@view inside_total[index_inside+1:end, 2])
    outside = sum(@view outside_total[index_outside+1:end, 2])
    tip, inside, outside
end

function extract_contact_forces(fcᵢ, grid, pile)
    inside = Dict{Float64, Float64}()
    outside = Dict{Float64, Float64}()
    for I in eachindex(fcᵢ) # walk from lower height
        y = grid[I][2] - tip_height(pile)
        fcy = -2π * fcᵢ[I][2]
        iszero(fcy) && continue
        centerline1 = Line(((getline(pile, 1) + reverse(getline(pile, 5))) / 2)...)
        centerline2 = Line(((getline(pile, 2) + reverse(getline(pile, 4))) / 2)...)
        isinside = GeometricObjects.ray_casting_to_right(centerline1, grid[I]) ||
        GeometricObjects.ray_casting_to_right(centerline2, grid[I])
        output = isinside ? inside : outside
        output[y] = get(output, y, 0.0) + fcy
    end
    function dict2array(dict)
        data′ = collect(dict)
        sort!(data′, by = x -> x[1])
        data = Array{Float64}(undef, length(data′), 2)
        @inbounds for i in 1:length(data′)
            data[i,1] = data′[i][1]
            data[i,2] = data′[i][2]
        end
        data
    end
    dict2array(inside), dict2array(outside)
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
