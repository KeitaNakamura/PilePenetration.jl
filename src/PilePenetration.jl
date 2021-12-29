module PilePenetration

using PoingrSimulator
using PoingrSimulator.Poingr
using PoingrSimulator.GeometricObjects
using PoingrSimulator.GeometricObjects: getline

using TOML

using Base: @_propagate_inbounds_meta, @_inline_meta

include("postprocess.jl")

function main_simulation()::Cint
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

function main_simulation(inputtoml_file::AbstractString)
    proj_dir = splitdir(inputtoml_file)[1]
    input = TOML.parsefile(inputtoml_file)

    General  = input["General"]
    Advanced = input["Advanced"]
    Pile     = input["Pile"]

    inputtoml = Dict{String, Any}(
        "General" => Dict{String, Any}(
            "simulation"        => "PenetrateIntoGround",
            "coordinate_system" => "axisymmetric",           # plane_strain or axisymmetric
            "domain"            => General["domain"],        # [[xmin, xmax], [ymin, ymax]] (m)
            "grid_space"        => General["grid_space"],    # (m)
            "gravity"           => General["gravity"],       # (m/s2)
            "total_time"        => General["total_time"],    # (sec)
            "show_progress"     => General["show_progress"],
        ),
        "Output" => Dict{String, Any}(
            "interval"       => General["output_interval"],   # (sec)
            "folder_name"    => joinpath(proj_dir, General["output_folder_name"]),
            "history"        => false,
            "serialize"      => true,
            "paraview"       => true,
            "paraview_grid"  => false,
            "copy_inputfile" => false,
        ),
        "Advanced" => Dict{String, Any}(
            "CFL"                       => Advanced["CFL"],
            "contact_threshold_scale"   => Advanced["contact_threshold_scale"],
            "contact_penalty_parameter" => Advanced["contact_penalty_parameter"],
            "npoints_in_cell"           => 2, # in each dimension
        ),
    )

    inputtoml["SoilLayer"] = map(input["SoilLayer"]) do SoilLayer
        friction_with_pile_inside  = SoilLayer["friction_with_pile_inside"]
        friction_with_pile_outside = SoilLayer["friction_with_pile_outside"]
        Dict{String, Any}(
            "thickness"               => SoilLayer["thickness"],                        # (m)
            "density"                 => SoilLayer["unit_weight"] / General["gravity"], # (kg/m3)
            "poissons_ratio"          => SoilLayer["poissons_ratio"],
            "youngs_modulus"          => SoilLayer["youngs_modulus"],                   # (N/m2)
            "cohesion"                => SoilLayer["cohesion"],                         # (N/m2)
            "friction_angle"          => SoilLayer["friction_angle"],                   # (degree)
            "dilatancy_angle"         => SoilLayer["dilatancy_angle"],                  # (degree)
            "tension_cutoff"          => SoilLayer["tension_cutoff"],
            "friction_with_rigidbody" => [friction_with_pile_inside, friction_with_pile_inside,
                                          friction_with_pile_outside, friction_with_pile_outside,
                                          friction_with_pile_outside, friction_with_pile_outside],
        )
    end

    inputtoml["RigidBody"] = begin
        radius_head    = Pile["diameter_head"] / 2 # (m)
        radius_tip     = Pile["diameter_tip"] / 2  # (m)
        pile_length    = Pile["pile_length"]       # (m)
        tapered_length = Pile["tapered_length"]    # (m)
        thickness      = Pile["thickness"]         # (m) NOTE: This must be at least `2 * grid_space`
        Dict{String, Any}(
            "velocity"    => Pile["velocity"],  # (m/s)
            "coordinates" => [[radius_head, pile_length], [radius_head, tapered_length], [radius_tip, 0.0],
                              [radius_tip+thickness, 0.0], [radius_head+thickness, tapered_length], [radius_head+thickness, pile_length]]
        )
    end

    inputtoml["Pile"] = input["Pile"]

    output_dir = inputtoml["Output"]["folder_name"]
    mkpath(output_dir)
    cp(inputtoml_file, joinpath(output_dir, "input.toml"); force = true)

    io = IOBuffer()
    TOML.print(io, inputtoml)
    PoingrSimulator.main(proj_dir, String(take!(io)), include(joinpath(@__DIR__, "injection.jl")))
end

end # module
