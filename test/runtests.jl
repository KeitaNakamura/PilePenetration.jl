using Test
using PilePenetration

using TOML
using CSV
using Serialization

@testset "Short pile penetration" begin
    PilePenetration.main_simulation("input.toml"); println()
    output_dir = TOML.parsefile("input.toml")["General"]["output_folder_name"]

    # check results
    output = CSV.File("output.csv") # expected output
    history = CSV.File(joinpath(output_dir, "history.csv"))
    for name in propertynames(output)
        output_col = output[name]
        history_col = history[name]
        @test history_col ≈ output_col atol=1e-4
    end

    # check input.toml
    @test TOML.parsefile("input.toml") == TOML.parsefile(joinpath(output_dir, "input.toml"))

    # check output files
    nums = 1:11
    for (name, ext) in (("paraview/output", ".vtm"),
                        ("force_inside/force_inside_", ".csv"),
                        ("force_outside/force_outside_", ".csv"))
        @test all(i -> isfile(joinpath(output_dir, string(name, i, ext))), nums)
    end
    for i in nums
        data = deserialize(joinpath(output_dir, "serialize", "save$i"))
        @test keys(data) === (:pointstate, :grid, :pile, :t)
    end
end
