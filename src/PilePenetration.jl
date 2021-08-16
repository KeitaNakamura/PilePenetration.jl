module PilePenetration

using Poingr, GeometricObjects
using DelimitedFiles

using Base: @_propagate_inbounds_meta, @_inline_meta

showval(val, name::String, unit::String = "") = println(rpad(name * ":", 10, '　') * lpad(string(val), 10) * " (" * unit * ")")
macro showval(ex, name, unit)
    quote
        v = $(esc(ex))
        showval(v, $(esc(name)), $unit)
    end
end
macro showval(ex, name)
    quote
        v = $(esc(ex))
        showval(v, $(esc(name)))
    end
end

struct NodeState
    f::Vec{2, Float64}
    fc::Vec{2, Float64}
    d::Vec{2, Float64}
    w::Float64
    m::Float64
    v::Vec{2, Float64}
    v_pile::Vec{2, Float64}
    vᵣ::Vec{2, Float64}
end

struct PointState
    m::Float64
    V0::Float64
    h::Vec{2, Float64}
    x::Vec{2, Float64}
    v::Vec{2, Float64}
    b::Vec{2, Float64}
    σ::SymmetricSecondOrderTensor{3, Float64, 6}
    F::SecondOrderTensor{3, Float64, 9}
    ∇v::SecondOrderTensor{3, Float64, 9}
    C::Mat{2, 3, Float64, 6}
end

function main()
    coord_system = :axisymmetric

    res = 1 # 数値が高いほど解像度が大きい
    @showval ρ₀ = 1.6e3      "乾燥密度"          "kg/m³"
    @showval g = 9.81        "重力加速度"        "m/s²"
    @showval h = 20          "地盤高さ"          "m"
    @showval ϕ = 35          "内部摩擦角"        "˚"
    @showval ψ = 0           "ダイレタンシー角"  "˚"
    @showval ν = 0.333       "ポアソン比"
    @showval E = 1e6         "ヤング係数"        "Pa"
    @showval μ = 0.3         "摩擦係数"
    # @showval dx = 0.1/res   "メッシュの幅"      "m"
    @showval dx = 0.025/res   "メッシュの幅"      "m"
    # @showval dx = 0.02/res   "メッシュの幅"      "m"

    @showval thick = 2dx           "肉厚"            "m"
    @showval D_i = 1.5             "杭頭径（内径）"  "m"
    @showval d_i = 1.221           "先端径（内径）"  "m"
    @showval taper_length = 4.0    "テーパー長"      "m"
    @showval taper_angle = atan((D_i-d_i) / 2taper_length) |> rad2deg "テーパー角" "˚"

    @showval vy_pile = 1.0         "貫入速度"        "m/s"
    @showval total_time = 14.8     "貫入時間"        "s"

    @showval process_num = 1  (process_num == 1 ? "吸い出しなし"     :
                               process_num == 2 ? "貫入中に吸い出し" :
                               process_num == 3 ? "貫入後に吸い出し" : error())
    @showval excavate_height = 1*D_i "吸い出し高さ" "m"
    if process_num == 3
        penetration_time = total_time
        total_time *= 2
    end

    grid = Grid(NodeState, LinRange(0:dx:6.0), LinRange(0:dx:40.0))
    pointstate = generate_pointstate((x,y) -> y < h, PointState, grid, coord_system)
    space = MPSpace(WLS{1}(CubicBSpline{2}()), grid, pointstate.x)
    elastic = LinearElastic(E = E, ν = ν)
    model = DruckerPrager(elastic, :inscribed, c = 0, ϕ = ϕ, ψ = ψ)

    @. pointstate.m = ρ₀ * pointstate.V0
    @. pointstate.F = one(SecondOrderTensor{3,Float64})
    @. pointstate.b = Vec(0.0, -g)

    for p in 1:length(pointstate)
        y = pointstate.x[p][2]
        σ_y = -ρ₀ * g * (h - y)
        σ_x = σ_y * ν / (1 - ν)
        pointstate.σ[p] = (@Mat [σ_x 0.0 0.0
                                 0.0 σ_y 0.0
                                 0.0 0.0 σ_x]) |> symmetric
    end

    R_i = D_i / 2 # radius
    r_i = d_i / 2 # radius
    R_o = R_i + thick
    r_o = r_i + thick
    tip_height_0 = h + gridsteps(grid, 1)
    y_max = grid[end, end][2]
    pile = Polygon([Vec(R_i, y_max),
                    Vec(R_i, tip_height_0 + taper_length),
                    Vec(r_i, tip_height_0),
                    Vec(r_o, tip_height_0),
                    Vec(R_o, tip_height_0 + taper_length),
                    Vec(R_o, y_max)])
    v_pile = Vec(0.0, -vy_pile)

    pile_center_0 = centroid(pile)

    find_ground_pos(xₚ) = maximum(x -> x[2], filter(x -> x[1] < gridsteps(grid, 1), xₚ))
    ground_pos0 = find_ground_pos(pointstate.x)


    @showval length(pointstate) "粒子数"

    # Output files
    ## proj
    proj_dir = joinpath("pile_finemesh.tmp")
    mkpath(proj_dir)

    ## paraview
    paraview_file = joinpath(proj_dir, "out")
    paraview_collection(vtk_save, paraview_file)

    ## history
    csv_file = joinpath(proj_dir, "history.csv")
    open(csv_file, "w") do io
        writedlm(io, ["disp" "force" "disp_inside_pile" "tip_resistance" "inside_resistance" "outside_resistance"], ',')
    end

    ## forces
    mkpath(joinpath(proj_dir, "force_tip"))
    mkpath(joinpath(proj_dir, "force_inside"))
    mkpath(joinpath(proj_dir, "force_outside"))

    ## copy this file
    cp(@__FILE__, joinpath(proj_dir, "main.jl"), force = true)

    t = 0.0
    logger = Logger(0.0:0.05:total_time; progress = true)
    # logger = Logger(0.0:0.0025:0.01; progress = true)
    while !isfinised(logger, t)
        reinit!(space, grid, pointstate.x, exclude = x -> in(x, pile))

        # ρ_min = minimum(@. pointstate.m / (pointstate.V0 * det(pointstate.F)))
        # vc = soundspeed(model.elastic.K, model.elastic.G, ρ_min)
        # dt = max(0.5 * minimum(gridsteps(grid)) / vc, 1e-5)
        # @show dt = 0.4*minimum(gridsteps(grid)) / vc
        dt = 1.0e-4

        default_point_to_grid!(grid, pointstate, space, coord_system)
        @. grid.state.v += (grid.state.f / grid.state.m) * dt

        if any(space.nearsurface)
            point_to_grid!((grid.state.d, grid.state.vᵣ), space, space.nearsurface) do it, p, i
                @_inline_meta
                @_propagate_inbounds_meta
                N = it.N
                w = it.w
                vₚ = pointstate.v[p]
                dₚ = contact_distance(pile, pointstate.x[p], pointstate.h[p])
                d = N * dₚ
                vᵣ = w * (vₚ - v_pile)
                d, vᵣ
            end
            @. grid.state.vᵣ /= grid.state.w
            @. grid.state.fc = contact_force(grid.state.vᵣ, grid.state.d, grid.state.m, dt, μ)
            @. grid.state.v += (grid.state.fc / grid.state.m) * dt
        end

        for bd in eachboundary(grid)
            @. grid.state.v[bd.indices] = boundary_velocity(grid.state.v[bd.indices], bd.n)
        end

        default_grid_to_point!(pointstate, grid, space, dt, coord_system)
        @inbounds Threads.@threads for p in eachindex(pointstate)
            Fₚ = pointstate.F[p]
            if det(Fₚ) < 0
                pointstate.F[p] = zero(Fₚ)
            end
            pointstate.σ[p] = stress(model, pointstate.σ[p], pointstate.∇v[p], dt)
        end

        translate!(pile, v_pile * dt)
        update!(logger, t += dt)

        if islogpoint(logger)
            Poingr.reordering_pointstate!(pointstate, space)
            paraview_collection(paraview_file, append = true) do pvd
                vtk_multiblock(string(paraview_file, logindex(logger))) do vtm
                    vtk_points(vtm, pointstate.x) do vtk
                        ϵₚ = @dot_lazy symmetric(pointstate.F - $Ref(I))
                        vtk_point_data(vtk, pointstate.v, "velocity")
                        vtk_point_data(vtk, @dot_lazy(-mean.(pointstate.σ)), "mean stress")
                        vtk_point_data(vtk, @dot_lazy(deviatoric_stress.(pointstate.σ)), "deviatoric stress")
                        vtk_point_data(vtk, @dot_lazy(volumetric_strain.(ϵₚ)), "volumetric strain")
                        vtk_point_data(vtk, @dot_lazy(deviatoric_strain.(ϵₚ)), "deviatoric strain")
                        vtk_point_data(vtk, pointstate.σ, "stress")
                        vtk_point_data(vtk, ϵₚ, "strain")
                        vtk_point_data(vtk, @dot_lazy(pointstate.m / (pointstate.V0 * det(pointstate.F))), "density")
                    end
                    vtk_grid(vtm, pile)
                    # vtk_grid(vtm, grid) do vtk
                        # vtk_point_data(vtk, vec(vᵢ), "nodal velocity")
                        # vtk_point_data(vtk, vec(vᵢ_before_contact), "nodal velocity before contact")
                        # vtk_point_data(vtk, vec(fcᵢ), "nodal contact force")
                    # end
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

            open(io -> writedlm(io, tip, ','), joinpath(proj_dir, "force_tip", "force_tip_$(logindex(logger)).csv"), "w")
            open(io -> writedlm(io, inside, ','), joinpath(proj_dir, "force_inside", "force_inside_$(logindex(logger)).csv"), "w")
            open(io -> writedlm(io, outside, ','), joinpath(proj_dir, "force_outside", "force_outside_$(logindex(logger)).csv"), "w")

            if process_num == 2 || process_num == 3
                if process_num == 2
                    height = (tip_height_0 - vy_pile * t) + excavate_height
                elseif process_num == 3
                    bottom_target = (tip_height_0 - vy_pile * penetration_time) + excavate_height
                    bottom_current = (h - bottom_target) / (total_time-penetration_time) * (t - penetration_time)
                    height = h - bottom_current
                end
            end
        end
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

function stress(model, σₚ, ∇vₚ, dt)
    σ = update_stress(model, σₚ, symmetric(∇vₚ) * dt)
    σ = Poingr.jaumann_stress(σ, σₚ, ∇vₚ, dt)
    mean(σ) > 0 ? zero(σ) : σ
end

function contact_distance(poly::Polygon, x::Vec{dim, T}, l::Vec{dim, T}) where {dim, T}
    threshold = mean(l) / 2 * 3
    d = distance(poly, x, threshold)
    d === nothing && return zero(Vec{dim, T})
    norm_d = norm(d)
    n = d / norm_d
    (threshold - norm_d) * n
end

function contact_force(vᵣ::Vec, d::Vec, m::Real, dt::Real, μ::Real)
    ξ = 0.5
    iszero(d) && return zero(d)
    n = d / norm(d)
    vᵣ_n = (vᵣ ⋅ n) * n
    vᵣ_t = vᵣ - vᵣ_n
    f_n = (1-ξ) * m * (2d/dt) / dt
    # f_n = 1e7 * d
    # f_n = (1-ξ) * m * (d/dt + vᵣ_n) / dt
    aᵣ_t = vᵣ_t / dt
    f_t = m * aᵣ_t
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
plot((arr = readdlm("pile.tmp/history.csv", ',', skipstart=1); @show size(arr, 1); (arr[:,1], arr[:,2])), legend = false, xlims = (0,2), ylims = (0,60e3))
plot((arr = readdlm("pile.tmp/history.csv", ',', skipstart=1); @show size(arr, 1); (arr[:,1], arr[:,[2,4,5,6]])), legend = false, xlims = (0,2), ylims = (0,60e3))
=#

end # module
