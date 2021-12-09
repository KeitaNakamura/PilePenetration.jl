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
