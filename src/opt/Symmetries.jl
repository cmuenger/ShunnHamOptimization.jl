
#generators = [(3), (2,1), (1,1,1)]
#generators = [(4), (3,1), (2,2), (2,1,1), (1,1,1,1)]
#generators = [(5), (4,1), (3,2), (3,1,1), (2,2,1), (2,1,1,1), (1,1,1,1,1)]
#generators = [(6), (5,1), (4,2), (3,3), (4,1,1), (3,2,1), (2,2,2), (3,1,1,1), (2,2,1,1), (2,1,1,1,1,1), (1,1,1,1,1,1)]



#build parametric description of symetric group in barycentric coordinates
#n is tuple that describes the repetition of a coordinates
#ex. n=(2,1,1) => x=(η1,η1,η2,1-2η1-η2)
function build_lamdba(n)
    x = Array{Sym,1}()
    i = 0;
    for j in 1:length(n)-1
        for k in 1:n[j]
            push!(x,Sym("η$i"))
        end
        i+=1
    end
    last = 1//n[end]*(1-sum(x))
    for k in 1:n[end]
        push!(x,last)
    end
    return (x,i)
end

#build all parametric description of symmtric configuration describe by the tuples in generator
function build_all_groups(generator)
    λ = []
    n = []
    for g in generator
        (x,i)=build_lamdba(g)
        push!(λ,x)
        push!(n,i)
    end
    return (λ,n)
end

#Compute size of all symmtric configuration for given generators and dimension
function size_all_groups(generator,d)
    size = [factorial(d+1) for i in 1:length(generator)]

    for (j,n) in enumerate(generator)
        f = 1;
        for i in n
            f *= factorial(i)
        end
        size[j] /= f
    end
    return size
end


#Compute points based on sphere-closed-packed arrangment
function build_latticepoints(n,d)
    if n==0
        return [Sym[1//(d+1) for i in 1:d+1] ]
    end
    λ = []
    for i in 0:n
        for m in multiexponents(d,i)
            x = Array{Sym}(undef,d+1)
            x[1:end-1] = 1//n*m
            x[end] = (1-sum(x[1:end-1]))
            push!(λ,x)
        end
    end
    return λ
end



function get_group_for_lattice_points(N,D,groups)
    if N==0
        return []
    end
    id = Array{Int64,1}()
    for i in 0:N
        for m in multiexponents(D,i)
            x = Array{Sym}(undef,D+1)
            x[1:end-1] = 1//N*m
            x[end] = (1-sum(x[1:end-1]))
            push!(id,get_group(groups,x))
        end
    end

    return id
end


function extract_init(pa,groups,n,generator)

    init = []
    while length(pa) > 0
        id = pa[1][1]
        λ = groups[id]
        g = generator[id]
        nvar = n[id]

        vars = []
        offset = 1
        for j in 1:length(g)-1
            push!(vars,pa[1][2][offset])
            offset += g[j]
        end
        #println(vars)

        push!(init,(id,vars))
        #println(pa)
        η_0 = [Sym("η$(l)") for l in 0:length(vars)-1 ]
        subn = [(η_0[l],vars[l]) for l in 1:length(vars)]

        λn = multiset_permutations(λ,length(λ))

        for λ_i in λn
            y =  λ_i.subs(subn)
            deleteat!(pa, findall(x->norm(x[2]-y)<1e-10,pa))
        end
    end 

    res = []
    for b in sort(init)
        for x in b[2]
            push!(res,x)
        end 
    end

    
    return res
end
function extract_init_from_lattticepoints(N,D,groups,generator,n,p)
    if N==0
        return []
    end
    η = Array{Float64}(undef,sum(n))
    idx = 1;
    for i in 0:N
        for m in multiexponents(D,i)
            x = Array{Sym}(undef,D+1)
            x[1:end-1] = 1//N*m
            x[end] = (1-sum(x[1:end-1]))
            println(groups)
            println(x)
            id = get_group_unique(groups,x)
            println("Group idx: ",id)

            if id > 1
                g=generator[id];
                offset = sum(n[1:id-1])
                offset2 = 1
                for j in 1:length(g)-1
                    η[offset+j] = p[idx][offset2]
                    offset2 += g[j]
                end
            end

            idx+=1
        end
    end
    return η
end

#Find number of groups for given parametric description for a certain symmetry groups
function build_supergroups(groups,points)

    sg = zeros(Int64,length(groups))
    for p in points
        i = get_group(groups,p)
        sg[i] += 1
    end

    return sg
end

#Find symmerty group for a parametric description
function get_group(groups,p)

   
    for (i,g) in enumerate(groups)
        group = multiset_permutations(g,length(g))
        for gp in group
            #println(gp)
            #println(g)
            if linsolve(tuple(gp-p...),(η0,η1,η2,η3,η4)) != ∅
               # println(gp)
                #println(p)
                #println(linsolve(tuple(gp-p...),(η0,η1,η2,η3,η4)) )
                return i
            end
        end
    end
    return 0
    
end

function get_group_unique(groups,p)

   
    for (i,g) in enumerate(groups)
    
            if linsolve(tuple(g-p...),(η0,η1,η2,η3,η4)) != ∅
               # println(gp)
                #println(p)
                #println(linsolve(tuple(gp-p...),(η0,η1,η2,η3,η4)) )
                return i
            end
        
    end
    
    println("Should not happen!!!")
end

function build_paraquad(supergroups,pgroups,n_p_group,s::Simplex)

    x = []
    w = []
    α = []
    γ = []

    #current index of α
    i = 0
    #current index of γ
    j = 0
    for (idx,sg) in enumerate(supergroups)
        for k in 1:sg
            n = n_p_group[idx]
      

            η_0 = [Sym("η$(l)") for l in 0:n-1 ]
            α_i = [Sym("α$(l+i+1)") for l in 0:n-1 ]

            if n == 0
                α_i  = []
            end
           
            γ_i = Sym("γ$(j+1)")

            suba = [(η_0[l],α_i[l]) for l in 1:n]


            
            λ = multiset_permutations(pgroups[idx],length(pgroups[idx]))

            push!(α, α_i)
            push!(γ, γ_i)

            for λ_i in λ
                push!(x,expand.(computeCartCoord(λ_i.subs(suba),s)))
                push!(w,γ_i)
            end

            i += n
            j += 1

        end
    end
    n = (i,j)
    p = Parameter(α,γ,n)
    return ParaQuad(x,w,p)
end


#Compute points and weights in barycentric coordinates
#Given super groups configuration, parametric description, count of free variables 
#values for parameter (η,γ) 
function print_pts_wts(supergroups, pgroups, n_g, η, γ)
    #current index of η
    i = 0
    j = 1
    for (idx,sg) in enumerate(supergroups)
        for k in 1:sg
            n = n_g[idx]
         
            η_0 = [Sym("η$(l)") for l in 0:n-1 ]
          
            suba = [(η_0[l],η[l+i]) for l in 1:n]
   
            λ = multiset_permutations(pgroups[idx],length(pgroups[idx]))
   
            for λ_i in λ
                println(λ_i.subs(suba)," ", γ[j])
            end
   
            i += n
            j += 1
           
       end
    end
end

function solve_para(groups,p)
    for (i,g) in enumerate(groups)
        group = multiset_permutations(g,length(g))
        for gp in group
    
            if linsolve(tuple(gp-p...),(η0,η1,η2,η3,η4)) != ∅
                #println(p)
                #println(gp)
                #println(linsolve(tuple(gp-p...),(η0,η1,η2,η3,η4)) )
                return;
            end
        end
    end
    println("Should not happen!!!")
end

#Compute innitial 
function compute_init(groups,points)
    for p in points
        solve_para(groups,p)
    end
end
