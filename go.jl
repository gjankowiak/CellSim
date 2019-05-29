try
    using Revise
catch
    import Pkg
    Pkg.add("Revise")
    using Revise
end

function go(config::String="configs/test_nucleus.yaml")
    println("Using configuration: '$config'")
    if length(ARGS) == 0
        push!(ARGS, config)
    else
        ARGS[1] = config
    end

    include("run.jl")
end

