module fet
    using SpecialFunctions: lgamma
    include("Fisher2x2.jl")
    export FisherExact2x2Test
end # module fet