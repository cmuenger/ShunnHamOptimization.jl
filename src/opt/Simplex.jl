# D; Dimension
# N; Num Vertices
# T: Type
struct Simplex{N,D,T}
    vertices::SVector{N,SVector{D,T}}
    bounds::SVector{D,SVector{2,T}}
    volume::T
end


function RefSimplex2d()

    a = √Sym(3)//6
    
    v = SVector{3,SVector{2,Sym}}(
         [ -δ//2 -√Sym(3)*δ//6 ],
         [ δ//2 -√Sym(3)*δ//6 ],
         [ 0 √Sym(3)*δ//3 ])

    b = SVector{2,SVector{2,Sym}}(
         [-√Sym(3)/3*(√Sym(3)/3*δ-Sym("x2"))  √Sym(3)/3*(√Sym(3)/3*δ-Sym("x2")) ],
         [-δ*√Sym(3)/6 δ*√Sym(3)/3])

    vol =  δ^2*√Sym(3)//4

    return Simplex(v,b,vol)
end

function RefSimplex3d()

    v = SVector{4,SVector{3,Sym}}(
        δ*[ -1//2 -√Sym(3)//6 -√Sym(6)//12],
        δ*[  1//2 -√Sym(3)//6 -√Sym(6)//12],
        δ*[    0   √Sym(3)//3 -√Sym(6)//12],
        δ*[    0       0       √Sym(6)//4])

    b = SVector{3,SVector{2,Sym}}(
        [ -√Sym(3)/3*(√Sym(2)/2*(δ*√Sym(6)/4-Sym("x3"))-Sym("x2")) √Sym(3)/3*(√Sym(2)/2*(δ*√Sym(6)/4-Sym("x3"))-Sym("x2"))],
        [ -√Sym(2)/4*(δ*√Sym(6)/4-Sym("x3")) √Sym(2)/2*(δ*√Sym(6)/4-Sym("x3"))],
        [ -δ*√Sym(6)/12 δ*√Sym(6)/4 ])
       
    vol =  δ^3*√Sym(2)//12

    return Simplex(v,b,vol)
end

function RefSimplex4d()

    v = SVector{5,SVector{4,Sym}}(
        δ*[ -1//2 -√Sym(3)/6 -√Sym(6)/12 -√Sym(10)/20],
        δ*[  1//2 -√Sym(3)/6 -√Sym(6)/12 -√Sym(10)/20],
        δ*[    0   √Sym(3)/3 -√Sym(6)/12 -√Sym(10)/20],
        δ*[    0        0     √Sym(6)/4  -√Sym(10)/20],
        δ*[    0        0          0      √Sym(10)/5])
        

    b = SVector{4,SVector{2,Sym}}(
        [ -√Sym(3)/3*(√Sym(2)/2*(√Sym(15)/5*(δ*√Sym(10)/5-Sym("x4"))-Sym("x3"))-Sym("x2")) √Sym(3)/3*(√Sym(2)/2*(√Sym(15)/5*(δ*√Sym(10)/5-Sym("x4"))-Sym("x3"))-Sym("x2"))],
        [ -√Sym(2)/4*(√Sym(15)/5*(δ*√Sym(10)/5-Sym("x4"))-Sym("x3")) √Sym(2)/2*(√Sym(15)/5*(δ*√Sym(10)/5-Sym("x4"))-Sym("x3"))],
        [ -√Sym(15)/15*(δ*√Sym(10)/5-Sym("x4")) √Sym(15)/5*(δ*√Sym(10)/5-Sym("x4")) ],
        [ -δ*√Sym(10)/20 δ*√Sym(10)/5 ])
       
    vol =  δ^4*√Sym(5)/96

    return Simplex(v,b,vol)
end

function RefSimplex5d()

    v = SVector{6,SVector{5,Sym}}(
        δ*[ -1//2 -√Sym(3)/6 -√Sym(6)/12 -√Sym(10)/20 -√Sym(15)/30 ],
        δ*[  1//2 -√Sym(3)/6 -√Sym(6)/12 -√Sym(10)/20 -√Sym(15)/30 ],
        δ*[    0   √Sym(3)/3 -√Sym(6)/12 -√Sym(10)/20 -√Sym(15)/30 ],
        δ*[    0        0     √Sym(6)/4  -√Sym(10)/20 -√Sym(15)/30 ],
        δ*[    0        0          0      √Sym(10)/5  -√Sym(15)/30 ],
        δ*[    0        0          0            0      √Sym(15)/6  ])

        

    b = SVector{5,SVector{2,Sym}}(
        [ -√Sym(3)/3*(√Sym(2)/2*(√Sym(15)/5*(√Sym(6)/3*(δ*√Sym(15)/6-Sym("x5"))-Sym("x4"))-Sym("x3"))-Sym("x2"))  √Sym(3)/3*(√Sym(2)/2*(√Sym(15)/5*(√Sym(6)/3*(δ*√Sym(15)/6-Sym("x5"))-Sym("x4"))-Sym("x3"))-Sym("x2")) ],
        [ -√Sym(2)/4*(√Sym(15)/5*(√Sym(6)/3*(δ*√Sym(15)/6-Sym("x5"))-Sym("x4"))-Sym("x3"))  √Sym(2)/2*(√Sym(15)/5*(√Sym(6)/3*(δ*√Sym(15)/6-Sym("x5"))-Sym("x4"))-Sym("x3")) ],
        [ -√Sym(15)/15*(√Sym(6)/3*(δ*√Sym(15)/6-Sym("x5"))-Sym("x4")) √Sym(15)/5*(√Sym(6)/3*(δ*√Sym(15)/6-Sym("x5"))-Sym("x4")) ],
        [ -√Sym(6)/12*(δ*√Sym(15)/6-Sym("x5")) √Sym(6)/3*(δ*√Sym(15)/6-Sym("x5")) ],
        [ -δ*√Sym(15)/30 δ*√Sym(15)/6 ])
       
    vol =  δ^5*√Sym(3)/480

    return Simplex(v,b,vol)
end

function integrate_f_simplex(f,s)

    integrand = f
    for (i,b) in enumerate(s.bounds)
        x = Sym("x$i")
        integrand =  integrate(integrand,(x,b[1],b[2]))
    end

    return integrand
end

function getBaryTransform(s)
    n= length(s.vertices)
    A = Array{Sym,2}(undef,n,n)

    for i in 1:n
         A[1:n-1,i] = s.vertices[i]
         A[n,i] = 1
    end
   
    SMatrix{n,n,Sym}(A)
end


function computeCartCoord(λ,s)

    sum(λ[i].*s.vertices[i] for i in 1:length(s.vertices))
end

function computeBaryCoord(x,s)

    A = getBaryTransform(s)
    An = A.evalf(subs=Dict(δ=>1.0))

    bn = [x.evalf(subs=Dict(δ=>1.0));1]
    λ = An\bn
    return λ
end