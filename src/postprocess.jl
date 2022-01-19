using DelimitedFiles
using PoingrSimulator: Input

function main_postprocess()::Cint
    inputtoml = only(ARGS)
    try
        main_postprocess(inputtoml)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function main_postprocess(inputtoml::AbstractString)
    proj_dir = splitdir(inputtoml)[1]
    INPUT = PoingrSimulator.parse_inputfile(inputtoml)
    postprocess_dir = joinpath(proj_dir, INPUT.General.output_folder_name)
    mkpath(postprocess_dir)
    cp(inputtoml, joinpath(postprocess_dir, "postprocess.toml"); force = true)
    main_postprocess(proj_dir, INPUT)
end

function main_postprocess(proj_dir::AbstractString, INPUT::Input{:Root})
    data = read_snapshots(joinpath(proj_dir, "snapshots.jld2"))
    postprocess_dir = joinpath(proj_dir, INPUT.General.output_folder_name)
    for name in keys(INPUT)
        if startswith(string(name), "Output")
            input = INPUT[name]
            @eval $(Symbol(lowercase(string(name))))($postprocess_dir, $input, $data)
        end
    end
end

function read_snapshots(path::AbstractString)
    data = jldopen(path, "r")
    map(i -> data[i], keys(data))
end

function outputhistory(postprocess_dir::AbstractString, INPUT::Input{:OutputHistory}, data)
    file = joinpath(postprocess_dir, "history.csv")

    ground_height_0 = INPUT.ground_height
    pile_center_0 = centroid(first(data).rigidbody0)

    outputhistory_head(file)
    for d in data
        outputhistory_append(
            file,
            d.grid,
            d.pointstate,
            d.rigidbody,
            INPUT.tip_height,
            INPUT.tapered_height,
            ground_height_0,
            pile_center_0,
        )
    end
end

function outputhistory_head(file::AbstractString)
    open(file, "w") do io
        write(io, join([
            "disp",
            "force",
            "disp_inside_pile",
            "tip",
            "inside",
            "outside",
            "taper",
            "straight",
            "tip_design",
            "inside_design",
            "outside_design",
            "taper_design",
            "straight_design",
        ], ",") * "\n")
    end
end

function outputhistory_append(
        file::AbstractString,
        grid,
        pointstate,
        pile::GeometricObject,
        tip_height,
        tapered_height,
        ground_height_0,
        pile_center_0,
    )
    inside_total, outside_total = extract_contact_forces(grid.state.fc, grid, pile[])
    open(file, "a") do io
        disp = norm(centroid(pile) - pile_center_0)
        force = -sum(grid.state.fc)[2] * 2π
        disp_inside_pile = -(find_ground_pos(pointstate.x, gridsteps(grid, 1)) - ground_height_0)
        tip, inside, outside = divide_force_into_tip_inside_outside(gridsteps(grid, 2), inside_total, outside_total)
        tip_design, inside_design, outside_design = divide_force_into_tip_inside_outside(tip_height, inside_total, outside_total)
        tip_taper, _, _ = divide_force_into_tip_inside_outside(tapered_height, inside_total, outside_total)
        taper = tip_taper - tip
        straight = force - tip_taper
        taper_design = tip_taper - tip_design
        straight_design = force - tip_taper
        write(io, join([
            disp,
            force,
            disp_inside_pile,
            tip,
            inside,
            outside,
            taper,
            straight,
            tip_design,
            inside_design,
            outside_design,
            taper_design,
            straight_design,
        ], ",") * "\n")
    end
end

function divide_force_into_tip_inside_outside(height_thresh, inside_total, outside_total)
    index_inside = searchsortedfirst(view(inside_total, :, 1), height_thresh) - 1
    index_outside = searchsortedfirst(view(outside_total, :, 1), height_thresh) - 1
    tip = sum(@view inside_total[begin:index_inside, 2]) + sum(@view outside_total[begin:index_outside, 2])
    inside = sum(@view inside_total[index_inside+1:end, 2])
    outside = sum(@view outside_total[index_outside+1:end, 2])
    tip, inside, outside
end

function extract_contact_forces(fcᵢ, grid, pile::Polygon)
    line_coordinates(pile, i) = coordinates(getline(pile, i))
    inside = Dict{Float64, Float64}()
    outside = Dict{Float64, Float64}()
    for I in eachindex(fcᵢ) # walk from lower height
        fcy = -2π * fcᵢ[I][2]
        iszero(fcy) && continue
        if grid[I][2] > pile[3][2] # straight and tapered part
            centerline1 = Line(((line_coordinates(pile, 1) + reverse(line_coordinates(pile, 5))) / 2)...)
            centerline2 = Line(((line_coordinates(pile, 2) + reverse(line_coordinates(pile, 4))) / 2)...)
            isinside = GeometricObjects.ray_casting_to_right(centerline1, grid[I]) ||
                       GeometricObjects.ray_casting_to_right(centerline2, grid[I])
        else # below the tip
            center = (pile[3][1] + pile[4][1]) / 2
            isinside = grid[I][1] < center
        end
        output = isinside ? inside : outside
        y = grid[I][2] - tip_height(pile)
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

tip_height(pile) = pile[3][2]
find_ground_pos(xₚ, dx) = maximum(x -> x[2], filter(x -> x[1] < dx, xₚ))
