struct Parameter
    α
    γ
    n
end

Parameter(n) = Parameter([Sym("α$i") for i in 1:n[1]],  [Sym("γ$i") for i in 1:n[2]], n)

function eval_parameter(p::Parameter, ϕ)

    α = (Sym("α$j")=>ϕ[j] for j in 1:p.n[1])
    γ = (Sym("γ$j")=>ϕ[j+p.n[1]] for j in 1:p.n[2])

    return (α,γ)
end

function get_alphas(p::Parameter)
    [Sym("α$i") for i in 1:p.n[1]]
end

function get_gammas(p::Parameter)
    [Sym("γ$i") for i in 1:p.n[2]]
end

function get_phi(p::Parameter)
    vcat(get_alphas(p),get_gammas(p))
end


struct ParaQuad
    x
    w
    p::Parameter
end 


function compute_quad(s,pQ,ϕ)

    A = getBaryTransform(s)
  
    An = A.evalf(subs=Dict(δ=>1.0))

    α,γ = eval_parameter(pQ.p,ϕ)

    for i in 1:length(pQ.x)

        x = pQ.x[i]
        w = pQ.w[i]

        b = [x.evalf(subs=Dict(δ=>1.0,α...));1]

        a1 = An\b 
        a2 = w(γ...)

        println(a1," ",a2)
    end
end

#=
#= 2D Rules =#
function parametrized_points_and_weights_3_2d(s)
  
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]

    x = [α1*e1, α1*e2, α1*e3] 
    w = [γ1 γ1 γ1 ]

    α = [α1]
    γ = [γ1]
    n=(1,1)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end

function parametrized_points_and_weights_6_2d(s)
 
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]

    x = [α1*e1, α1*e2, α1*e3,
         α2*(e1+e2), α2*(e1+e3), α2*(e2+e3)] 
    w = [γ1 γ1 γ1 γ2 γ2 γ2 ]

    α = [α1 α2]
    γ = [γ1 γ2]
    n=(2,2)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end

function parametrized_points_and_weights_10_2d(s)
    
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]

    x = [α1*e1, α1*e2, α1*e3,
         α2*e1+α3*e2, α2*e2+α3*e1,
         α2*e1+α3*e3, α2*e3+α3*e1,
         α2*e2+α3*e3, α2*e3+α3*e2,
         (e1+e2+e3)/3]
    w = [γ1 γ1 γ1 γ2 γ2 γ2 γ2 γ2 γ2 γ3 ]

    α = [α1, [α2 α3], 1]
    γ = [γ1 γ2 γ3]
    n=(3,3)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end

function parametrized_points_and_weights_15_2d(s)
  
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]

    x = [α1*e1, α1*e2, α1*e3,
         α2*(e1+e2), α2*(e1+e3), α2*(e2+e3),
         α3*e1+α4*e2, α3*e2+α4*e1,
         α3*e1+α4*e3, α3*e3+α4*e1,
         α3*e2+α4*e3, α3*e3+α4*e2,
         α5*e1, α5*e2, α5*e3]
    w = [γ1 γ1 γ1 γ2 γ2 γ2 γ3 γ3 γ3  γ3 γ3 γ3 γ4 γ4 γ4 ]

    α = [α1, α2, [α3 α4], α5]
    γ = [γ1 γ2 γ3 γ4]
    n=(5,4)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end
function parametrized_points_and_weights_21_2d(s)
  
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]

    x = [α1*e1, α1*e2, α1*e3,
         α2*(e1+e2), α2*(e1+e3), α2*(e2+e3),
         α3*e1+α4*e2, α3*e2+α4*e1,
         α3*e1+α4*e3, α3*e3+α4*e1,
         α3*e2+α4*e3, α3*e3+α4*e2,
         α5*e1+α6*e2, α5*e2+α6*e1,
         α5*e1+α6*e3, α5*e3+α6*e1,
         α5*e2+α6*e3, α5*e3+α6*e2,
         α7*e1, α7*e2, α7*e3]
    w = [γ1 γ1 γ1 γ2 γ2 γ2 γ3 γ3 γ3 γ3 γ3 γ3 γ4 γ4 γ4 γ4 γ4 γ4 γ5 γ5 γ5 ]

    α = [α1, α2, [α3 α4], [α5 α6], α7]
    γ = [γ1 γ2 γ3 γ4 γ5]
    n=(7,5)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end


#= 3D Rules =#
function parametrized_points_and_weights_4_3d(s)
   
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]
    e4 = s.vertices[4]

    x = [α1*e1, α1*e2, α1*e3, α1*e4 ] 
    w = [γ1 γ1 γ1 γ1]

    α = [α1]
    γ = [γ1]
    n=(1,1)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end

function parametrized_points_and_weights_10_3d(s)
    
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]
    e4 = s.vertices[4]

    x = [α1*e1, α1*e2, α1*e3, α1*e4,
         α2*(e1+e2), α2*(e1+e3), α2*(e1+e4), α2*(e2+e3),  α2*(e2+e4), α2*(e3+e4)] 
    w = [γ1 γ1 γ1 γ1 γ2 γ2 γ2 γ2 γ2 γ2]

    α = [α1 α2]
    γ = [γ1 γ2]
    n=(2,2)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end

function parametrized_points_and_weights_20_3d(s)
   
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]
    e4 = s.vertices[4]

    x = [α1*e1, α1*e2, α1*e3, α1*e4,
         α2*e1+α3*e2, α2*e2+α3*e1, 
         α2*e1+α3*e3, α2*e3+α3*e1,
         α2*e1+α3*e4, α2*e4+α3*e1, 
         α2*e2+α3*e3, α2*e3+α3*e2, 
         α2*e2+α3*e4, α2*e4+α3*e2, 
         α2*e3+α3*e4, α2*e4+α3*e3, 
         α4*(e1+e2+e3), α4*(e1+e2+e4), α4*(e1+e3+e4), α4*(e2+e3+e4)] 
    w = [γ1 γ1 γ1 γ1 γ2 γ2 γ2 γ2 γ2 γ2 γ2 γ2 γ2 γ2 γ2 γ2 γ3 γ3 γ3 γ3]

    α = [α1, [α2 α3], α4]
    γ = [γ1 γ2 γ3]
    n= (4,3)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end

function parametrized_points_and_weights_35_3d(s)
   
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]
    e4 = s.vertices[4]

    x = [α1*e1, α1*e2, α1*e3, α1*e4,
         α2*e1+α3*e2, α2*e2+α3*e1, 
         α2*e1+α3*e3, α2*e3+α3*e1,
         α2*e1+α3*e4, α2*e4+α3*e1, 
         α2*e2+α3*e3, α2*e3+α3*e2, 
         α2*e2+α3*e4, α2*e4+α3*e2, 
         α2*e3+α3*e4, α2*e4+α3*e3, 
         α4*(e1+e2),
         α4*(e1+e3),
         α4*(e1+e4),
         α4*(e2+e3),
         α4*(e2+e4),
         α4*(e3+e4),
         α5*e1+α6*(e2+e3),
         α5*e1+α6*(e2+e4),
         α5*e1+α6*(e3+e4),
         α5*e2+α6*(e1+e3),
         α5*e2+α6*(e1+e4),
         α5*e2+α6*(e3+e4),
         α5*e3+α6*(e2+e1),
         α5*e3+α6*(e2+e4),
         α5*e3+α6*(e1+e4),
         α5*e4+α6*(e2+e3),
         α5*e4+α6*(e2+e1),
         α5*e4+α6*(e3+e1),
         (e1+e2+e3+e4)//4] 
    w = [γ1, γ1, γ1, γ1,
         γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2,
         γ3, γ3, γ3, γ3, γ3, γ3,
         γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4,
         γ5]

    α = [α1, [α2 α3], α4,[α5 α6], 1]
    γ = [γ1 γ2 γ3 γ4 γ5]

    n = (6,5)

    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end

function parametrized_points_and_weights_56_3d()
    
    e1 = s.vertices[1]
    e2 = s.vertices[2]
    e3 = s.vertices[3]
    e4 = s.vertices[4]

    x = [α1*e1, α1*e2, α1*e3, α1*e4,
         α2*e1+α3*e2, α2*e2+α3*e1, 
         α2*e1+α3*e3, α2*e3+α3*e1,
         α2*e1+α3*e4, α2*e4+α3*e1, 
         α2*e2+α3*e3, α2*e3+α3*e2, 
         α2*e2+α3*e4, α2*e4+α3*e2, 
         α2*e3+α3*e4, α2*e4+α3*e3, 
         α4*e1+α5*e2, α4*e2+α5*e1, 
         α4*e1+α5*e3, α4*e3+α5*e1,
         α4*e1+α5*e4, α4*e4+α5*e1, 
         α4*e2+α5*e3, α4*e3+α5*e2, 
         α4*e2+α5*e4, α4*e4+α5*e2, 
         α4*e3+α5*e4, α4*e4+α5*e3, 
         α6*e1+α7*(e2+e3),
         α6*e1+α7*(e2+e4),
         α6*e1+α7*(e3+e4),
         α6*e2+α7*(e1+e3),
         α6*e2+α7*(e1+e4),
         α6*e2+α7*(e3+e4),
         α6*e3+α7*(e2+e1),
         α6*e3+α7*(e2+e4),
         α6*e3+α7*(e1+e4),
         α6*e4+α7*(e2+e3),
         α6*e4+α7*(e2+e1),
         α6*e4+α7*(e3+e1),
         α8*e1+α9*(e2+e3),
         α8*e1+α9*(e2+e4),
         α8*e1+α9*(e3+e4),
         α8*e2+α9*(e1+e3),
         α8*e2+α9*(e1+e4),
         α8*e2+α9*(e3+e4),
         α8*e3+α9*(e2+e1),
         α8*e3+α9*(e2+e4),
         α8*e3+α9*(e1+e4),
         α8*e4+α9*(e2+e3),
         α8*e4+α9*(e2+e1),
         α8*e4+α9*(e3+e1),
         α10*e1, α10*e2, α10*e3, α10*e4
    ]

    w = [γ1, γ1, γ1, γ1,
         γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2, γ2,
         γ3, γ3, γ3, γ3, γ3, γ3, γ3, γ3, γ3, γ3, γ3, γ3,
         γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4, γ4,
         γ5, γ5, γ5, γ5, γ5, γ5, γ5, γ5, γ5, γ5, γ5, γ5,
         γ6, γ6, γ6, γ6]

    α = [α1, [α2 α3], [α4 α5], [α6 α7], [α8 α9], α10]
    γ = [γ1 γ2 γ3 γ4 γ5 γ6]

    n = (10,6)
    
    p = Parameter( α, γ, n)
    return ParaQuad(x,w,p)
end
=#