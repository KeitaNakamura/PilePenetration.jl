module Injection

using PilePenetration
using PoingrSimulator.Poingr
using PoingrSimulator.GeometricObjects
using JLD2

const pile_center_0 = Ref(Vec(NaN, NaN))
const ground_height_0 = Ref(NaN)

function main_output_initialize(args)
    INPUT = args.INPUT

    if haskey(INPUT.General, :restart)
        grid, pointstate, pile = load(joinpath(INPUT.Output.original_directory, "snapshots.jld2"))["0"]
    else
        grid       = args.grid
        pointstate = args.pointstate
        pile       = args.rigidbody
    end

    pile_center_0[] = centroid(pile)
    ground_height_0[] = PilePenetration.find_ground_pos(pointstate.x, gridsteps(grid, 1))

    history_file = joinpath(INPUT.Output.directory, "history.csv")
    PilePenetration.outputhistory_head(history_file)
end

function main_output(args)
    grid         = args.grid
    pointstate   = args.pointstate
    pile         = args.rigidbody
    INPUT        = args.INPUT

    history_file = joinpath(INPUT.Output.directory, "history.csv")
    force_inside_directory = joinpath(INPUT.Output.directory, "force inside directory")
    force_outside_directory = joinpath(INPUT.Output.directory, "force outside directory")

    PilePenetration.outputhistory_append(
        history_file,
        grid,
        pointstate,
        pile,
        2 * pile[1][1], # 1D = 2R
        pile[2][2] - pile[3][2], # tapered length
        ground_height_0[],
        pile_center_0[],
    )

    if args.INPUT.Pile.vacuum
        vacuum_height = INPUT.Pile.vacuum_height
        pile_inside = Polygon(Vec(0.0, pile[1][2]), Vec(0.0, pile[3][2]), pile[3], pile[2], pile[1])
        inds = findall(xₚ -> (xₚ in pile_inside) && xₚ[2] > (PilePenetration.tip_height(pile) + vacuum_height), pointstate.x)
        deleteat!(pointstate, inds)
    end
end

end
