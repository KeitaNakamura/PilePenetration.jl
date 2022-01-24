using Test
using PilePenetration

using TOML
using CSV

using Serialization
using NaturalSort

@testset "Short pile penetration" begin
    PilePenetration.main_simulation("input.toml"); println()
    output_dir = TOML.parsefile("input.toml")["General"]["output_folder_name"]

    # check results
    output = CSV.File("output.csv") # expected output
    history = CSV.File(joinpath(output_dir, "history.csv"))
    for name in propertynames(output)
        output_col = output[name]
        history_col = history[name]
        @test history_col ≈ output_col atol=1e-3
    end

    # check input.toml
    @test TOML.parsefile("input.toml") == TOML.parsefile(joinpath(output_dir, "input.toml"))

    # check paraview files
    @test all(i -> isfile(joinpath(output_dir, string("paraview/output", i, ".vtm"))), 0:10)
    # check snapshots file
    root, _, files = only(walkdir(joinpath(output_dir, "snapshots")))
    count = 0
    for file in sort(files; lt = natural)
        @test file == "snapshot$count"
        data = deserialize(joinpath(root, file))
        @test data isa NamedTuple
        @test keys(data) == (:grid, :pointstate, :rigidbody, :rigidbody0, :t)
        count += 1
    end
end

@testset "Restart simulation" begin
    PilePenetration.main_simulation("input_restart.toml"); println()
    dict = TOML.parsefile("input_restart.toml")
    dir_ext = splitext(dict["General"]["output_folder_name"])
    output_dir = string(dir_ext[1], "_restarted_from_", dict["General"]["restart"], dir_ext[2])

    # check results
    output = CSV.File("output.csv") # expected output
    history = CSV.File(joinpath(output_dir, "history.csv"))
    for name in propertynames(output)
        output_last = output[name][end]
        history_last = history[name][end]
        @test history_last ≈ output_last atol=1e-3
    end
end

@testset "Post-processing" begin
    @testset "normal" begin
        output_dir = "output.tmp"
        inputfile = joinpath(output_dir, "postprocess.toml")
        cp("postprocess.toml", inputfile; force = true)
        PilePenetration.main_postprocess(inputfile)

        # check results
        output = CSV.File("output.csv") # expected output
        history = CSV.File(joinpath(output_dir, "postprocess.tmp", "history.csv"))
        for name in propertynames(output)
            output_col = output[name]
            history_col = history[name]
            @test history_col ≈ output_col atol=1e-3
        end
    end
    @testset "restart" begin
        output_dir = "output_restarted_from_5.tmp"
        inputfile = joinpath(output_dir, "postprocess.toml")
        cp("postprocess.toml", inputfile; force = true)
        PilePenetration.main_postprocess(inputfile)

        # check results
        output = CSV.File("output.csv") # expected output
        history = CSV.File(joinpath(output_dir, "postprocess.tmp", "history.csv"))
        for name in propertynames(output)
            output_last = output[name][end]
            history_last = history[name][end]
            @test history_last ≈ output_last atol=1e-3
        end
    end
end
