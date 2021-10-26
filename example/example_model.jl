using ShunnHamOptimization
using SymPy
using Plots
D = 2
N = 5

#Reference Simplex
s = RefSimplex2d()

#Build parametrized quadrature
(λ,n) = build_all_groups(generators2D)
n_g = size_all_groups(generators2D,D)
p = build_latticepoints(N,D)
sg = div.( build_supergroups(λ,p),n_g)

param2 = build_paraquad(sg,λ,n,s)

#Get initial points on lattice
p = build_latticepoints(N,D)

#2d: 3 | 3d: 6 | 4d: 10 | 5d: 15
#Scaling factor
a = N/(sqrt(3)+N)

pa = []
xa = []
for p_i in p 
    x = computeCartCoord(p_i,s)
    push!(pa,computeBaryCoord(a*x,s))
    push!(xa,a*x.evalf(subs=Dict(δ=>1.0)))
end

g_p = get_group_for_lattice_points(N,D,λ)
g_pa = collect(zip(g_p,pa))



#Plot scaled lattice points
#va = [v.evalf.(subs=Dict(δ=>1.0)) for v in s.vertices]
#scatter([v[1] for v in va], [v[2] for v in va])
#scatter!([a[1] for a in xa], [a[2] for a in xa])



init = extract_init(g_pa,λ,n,generators2D)


objective,constraints,order = getOptimizationEquations(D,s,param2)

println("======")
println("Order: ",order)
println("Objective:")
display(objective[1])
display(objective[2])
println("Contraints:")
for c in constraints
    display(c)
end


build_optmodel(objective,constraints,param2.p.n,init,length(p),"OptModelX.jl")

