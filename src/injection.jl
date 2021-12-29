module Injection

using PilePenetration
using PilePenetration.PoingrSimulator.Poingr
using PilePenetration.PoingrSimulator.GeometricObjects
using PilePenetration.DelimitedFiles

const pile_center_0 = Ref(Vec(NaN, NaN))
const ground_height_0 = Ref(NaN)

function main_output(args)
    output_index = args.output_index
    grid         = args.grid
    pointstate   = args.pointstate
    pile         = args.rigidbody
    INPUT        = args.INPUT

    history_file = joinpath(args.output_dir, "history.csv")
    force_inside_directory = joinpath(args.output_dir, "force inside directory")
    force_outside_directory = joinpath(args.output_dir, "force outside directory")

    if output_index == 0
        pile_center_0[] = centroid(pile)
        ground_height_0[] = PilePenetration.find_ground_pos(pointstate.x, gridsteps(grid, 1))
        PilePenetration.outputhistory_head(history_file)
    end

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
