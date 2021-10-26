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
