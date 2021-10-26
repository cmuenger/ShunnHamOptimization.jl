
using LinearAlgebra
using ShunnHamOptimization


#Load Model from file 
include("../models/2D/P15.jl")

n=1 #Number of optimizations restarts,
m=5 #Number of tested initial points 
model_j = P15_model_j


function find_opt(model)
    minf_g = Inf;
    minx_g = []

    for i in 1:100
    
        (minf,minx) = run_opt_sequential(model)
        r = test_constraints(model,minx)
        if norm(r) < 1e-10
            println("Min_x: $minx")
            println("constaints error $r")
        end

        if norm(r) < minf_g
            minf_g = norm(r)
            minx_g = minx
        end
    end

    return minx_g
end

minx = find_opt(model_j)

r = test_constraints(model_j, minx)
println("Min_x: $minx")
println("constaints error $r")

#Alternative version use SymPy expression directly (not recommended as it is slow)
#= TODO =#
