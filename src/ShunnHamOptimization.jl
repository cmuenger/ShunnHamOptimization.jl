module ShunnHamOptimization

using Combinatorics
using LinearAlgebra
using NLopt
using RowEchelon
using StaticArrays
using SymPy
using PyCall



export getOptimizationEquations

export generators2D,generators3D,generators4D,generators4D
const generators2D = [(3), (2,1), (1,1,1)]
const generators3D = [(4), (3,1), (2,2), (2,1,1), (1,1,1,1)]
const generators4D = [(5), (4,1), (3,2), (3,1,1), (2,2,1), (2,1,1,1), (1,1,1,1,1)]
const generators5D  = [(6), (5,1), (4,2), (3,3), (4,1,1), (3,2,1), (2,2,2), (3,1,1,1), (2,2,1,1), (2,1,1,1,1,1), (1,1,1,1,1,1)]

export build_all_groups,size_all_groups,build_latticepoints,build_supergroups,extract_init,print_pts_wts,build_paraquad, get_group_for_lattice_points

export Simplex
export RefSimplex2d,RefSimplex3d,RefSimplex4d,RefSimplex5d
export computeCartCoord,computeBaryCoord,compute_init

export Parameter,ParaQuad
export FunctionGradient
export OptModel,JuliaModel,numParameters,lower_bounds,upper_bounds,get_phi,build_optmodel
export build_model_julia,build_model_cxx

export print_pts_wts

export run_opt,run_opt_sequential,test_constraints

export δ,α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α11,α2,γ1,γ2,γ3,γ4,γ5,γ6,γ7

#=
const η0 = Sym("η0")
const η1 = Sym("η1")
const η2 = Sym("η2")
const η3 = Sym("η3")
const η4 = Sym("η4")

const ∅ = sympy.EmptySet()

#const δ = Sym("δ")

const α1 = Sym("α1")
const α2 = Sym("α2")
const α3 = Sym("α3")
const α4 = Sym("α4")
const α5 = Sym("α5")
const α6 = Sym("α6")
const α7 = Sym("α7")
const α8 = Sym("α8")
const α9 = Sym("α9")
const α10 = Sym("α10")
const α11 = Sym("α11")
const α12 = Sym("α12")


const γ1 = Sym("γ1")
const γ2 = Sym("γ2")
const γ3 = Sym("γ3")
const γ4 = Sym("γ4")
const γ5 = Sym("γ5")
const γ6 = Sym("γ6")
const γ7 = Sym("γ7")

#=
for op in [:δ]
    @eval begin
        const $op = Sym(C_NULL)
    end
    eval(Expr(:export, op))
end

macro init_constant(op, libnm)
    tup = (Base.Symbol("basic_const_$libnm"), libshunnhamoptimization)
    alloc_tup = (:basic_new_stack, libshunnhamoptimization)
    :(
        begin
            ccall($alloc_tup, Nothing, (Ref{Basic}, ), $op)
            ccall($tup, Nothing, (Ref{Basic}, ), $op)
            finalizer(basic_free, $op)
        end
    )
end

function init_constants()
    @init_constant δ Sym("δ")
end
=#

const η0 = Sym("η0")
const η1 = Sym("η1")
const η2 = Sym("η2")
const η3 = Sym("η3")
const η4 = Sym("η4")

const ∅ = sympy.EmptySet()

#const δ = Sym("δ")

const α1 = Sym("α1")
const α2 = Sym("α2")
const α3 = Sym("α3")
const α4 = Sym("α4")
const α5 = Sym("α5")
const α6 = Sym("α6")
const α7 = Sym("α7")
const α8 = Sym("α8")
const α9 = Sym("α9")
const α10 = Sym("α10")
const α11 = Sym("α11")
const α12 = Sym("α12")


const γ1 = Sym("γ1")
const γ2 = Sym("γ2")
const γ3 = Sym("γ3")
const γ4 = Sym("γ4")
const γ5 = Sym("γ5")
const γ6 = Sym("γ6")
const γ7 = Sym("γ7")
=#

pynull() = PyCall.PyNULL()
const δ = Sym(pynull())

const η0 = Sym(pynull())
const η1 = Sym(pynull())
const η2 = Sym(pynull())
const η3 = Sym(pynull())
const η4 = Sym(pynull())

const ∅ = Sym(pynull())

const α1 = Sym(pynull())
const α2 = Sym(pynull())
const α3 = Sym(pynull())
const α4 = Sym(pynull())
const α5 = Sym(pynull())
const α6 = Sym(pynull())
const α7 = Sym(pynull())
const α8 = Sym(pynull())
const α9 = Sym(pynull())
const α10 = Sym(pynull())
const α11 = Sym(pynull())
const α12 = Sym(pynull())


const γ1 = Sym(pynull())
const γ2 = Sym(pynull())
const γ3 = Sym(pynull())
const γ4 = Sym(pynull())
const γ5 = Sym(pynull())
const γ6 = Sym(pynull())
const γ7 = Sym(pynull())

function __init__()
    #init constant symbols
    copy!(δ.__pyobject__,  Sym("δ").__pyobject__)

    copy!(η0.__pyobject__,  Sym("η0").__pyobject__)
    copy!(η1.__pyobject__,  Sym("η1").__pyobject__)
    copy!(η2.__pyobject__,  Sym("η2").__pyobject__)
    copy!(η3.__pyobject__,  Sym("η3").__pyobject__)
    copy!(η4.__pyobject__,  Sym("η4").__pyobject__)

    copy!(∅.__pyobject__,  sympy.EmptySet().__pyobject__)

    copy!(α1.__pyobject__,  Sym("α1").__pyobject__)
    copy!(α2.__pyobject__,  Sym("α2").__pyobject__)
    copy!(α3.__pyobject__,  Sym("α3").__pyobject__)
    copy!(α4.__pyobject__,  Sym("α4").__pyobject__)
    copy!(α5.__pyobject__,  Sym("α5").__pyobject__)
    copy!(α6.__pyobject__,  Sym("α6").__pyobject__)
    copy!(α7.__pyobject__,  Sym("α7").__pyobject__)
    copy!(α8.__pyobject__,  Sym("α8").__pyobject__)
    copy!(α9.__pyobject__,  Sym("α9").__pyobject__)
    copy!(α10.__pyobject__,  Sym("α10").__pyobject__)
    copy!(α11.__pyobject__,  Sym("α11").__pyobject__)
    copy!(α12.__pyobject__,  Sym("α12").__pyobject__)

    copy!(γ1.__pyobject__,  Sym("γ1").__pyobject__)
    copy!(γ2.__pyobject__,  Sym("γ2").__pyobject__)
    copy!(γ3.__pyobject__,  Sym("γ3").__pyobject__)
    copy!(γ4.__pyobject__,  Sym("γ4").__pyobject__)
    copy!(γ5.__pyobject__,  Sym("γ5").__pyobject__)
    copy!(γ6.__pyobject__,  Sym("γ6").__pyobject__)
    copy!(γ7.__pyobject__,  Sym("γ7").__pyobject__)

end


include("opt/Simplex.jl")
include("opt/ParameterQaudrature.jl")
include("opt/Model.jl")
include("opt/Symmetries.jl")
include("opt/Constraints.jl")
include("opt/Optimization.jl")




end # module
