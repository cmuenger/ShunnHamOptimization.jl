using ShunnHamOptimization
using SymPy

# Load model
include("../models/2D/OptModel2D_28P.jl")

#generate julia code
name = "P28_model_j"
filename  = "P28.jl"
build_model_julia(P28_model_s,name,filename) 

#generate c++ code
name = "P28_model_c"
filename = "P28"
build_model_cxx(P28_model_s,name,filename) 
