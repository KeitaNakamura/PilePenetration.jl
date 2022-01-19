module Injection

using PilePenetration
using PoingrSimulator.Poingr
using PoingrSimulator.GeometricObjects
using JLD2

function main_output_initialize(args)
    INPUT = args.INPUT

    history_file = joinpath(INPUT.Output.directory, "history.csv")
    PilePenetration.outputhistory_head(history_file)
end

function main_output(args)
    grid         = args.grid
    pointstate   = args.pointstate
    pile         = args.rigidbody
    pile0        = args.rigidbody0
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
        sum(layer -> layer.thickness, INPUT.SoilLayer),
        centroid(pile0),
    )

    if args.INPUT.Pile.vacuum
        vacuum_height = INPUT.Pile.vacuum_height
        pile_inside = Polygon(Vec(0.0, pile[1][2]), Vec(0.0, pile[3][2]), pile[3], pile[2], pile[1])
        inds = findall(xₚ -> (xₚ in pile_inside) && xₚ[2] > (PilePenetration.tip_height(pile) + vacuum_height), pointstate.x)
        deleteat!(pointstate, inds)
    end
end

end
