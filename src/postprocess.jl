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

function outputhistory_append(file::AbstractString, grid, pointstate, pile, tip_height, tapered_height, ground_height_0, pile_center_0)
    inside_total, outside_total = extract_contact_forces(grid.state.fc, grid, pile)
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

function extract_contact_forces(fcᵢ, grid, pile)
    inside = Dict{Float64, Float64}()
    outside = Dict{Float64, Float64}()
    for I in eachindex(fcᵢ) # walk from lower height
        fcy = -2π * fcᵢ[I][2]
        iszero(fcy) && continue
        if grid[I][2] > pile[3][2] # straight and tapered part
            centerline1 = Line(((getline(pile, 1) + reverse(getline(pile, 5))) / 2)...)
            centerline2 = Line(((getline(pile, 2) + reverse(getline(pile, 4))) / 2)...)
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
