module PilePenetration

using Poingr
using GeometricObjects
using GeometricObjects: getline

using DelimitedFiles
using Serialization
using TOML
using Dates

using Base: @_propagate_inbounds_meta, @_inline_meta

function parseinput(dict::Dict)
    dict2namedtuple(x::Dict) = (; (Symbol(key) => value for (key, value) in x)...)
    list = map(collect(keys(dict))) do section
        content = dict[section]
        if content isa Dict
            Symbol(section) => dict2namedtuple(content)
        elseif content isa Vector
            Symbol(section) => map(dict2namedtuple, content)
        else
            error("unreachable")
        end
    end
    (; list...)
end
function parseinput(inputtoml::AbstractString)
    parseinput(TOML.parsefile(inputtoml))
end

include("simulation.jl")
include("postprocess.jl")

end # module
