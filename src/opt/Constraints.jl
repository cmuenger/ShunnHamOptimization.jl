
#Build symbolic partial derivative from multi index 'idx'
function partial_derivative(idx)
  
    s = String("")
   
    for (k,e) in enumerate(idx)
               
        s *= "x$k"^e
    end
           
    Sym("f_"*s)
end

#Build symbolic d-dimesional Taylor expansion around origin up to order n
function multitaylor(n,d)
    sum = 0

    for i in 0:n
       sum+= nth_multitaylor(i,d)
    end
    sum
end

#Build symbolic expression of the n-th order term of the d-dimesnional Taylor expansion
function nth_multitaylor(n,d)
    x =[]
    for i in 1:d
        push!(x,Sym("x$i"))
    end
   
    sum = 0

    multiindex = collect(multiexponents(d,n))
       
    for idx in multiindex
        prod = 1
            
        for (k,e) in enumerate(idx)
            prod *= (x[k]^e) 
        end
           
       
        f = partial_derivative(idx)
            
        sum+=1//factorial(n)*multinomial(idx...)*f*prod
    end

    sum
end

#Build symbolic expression for n-th degree polynom
function nth_deegree_poly(n,d)
    x =[]
    for i in 1:d
        push!(x,Sym("x$i"))
    end
   
    sum = 0

    exponents = collect(multiexponents(d,n))
       
    for expn in exponents
        prod = 1
            
        for (k,e) in enumerate(expn)
            prod *= (x[k]^e) 
        end
           
         #Constant factor
        f = partial_derivative(expn)
            
        sum+=f*prod
    end

    sum
end


function f_sub(d,x, w,f)

    res = 0.0
    for i in 1:length(x)
        subx = (Sym("x$j")=>x[i][j] for j in 1:d)

        res += w[i]*f(subx...)
    end
    res
end

#Compute error up to given order
function error(d,x,w,f,s)
    err = 0
    
    V = s.volume
    err = expand(integrate_f_simplex(f,s)/V-f_sub(d,x,w,f))
   
    err 
end

#Get reduced row echelon of the coefficient matrix and and vector u
function getCoeffMatrix(Cu,n,p)

    u = [Sym(1)]
    α = []
    for i in 1:length(p.α)
        if p.α[i] == 0
            push!(u,p.γ[i])
        else 
            for k in 0:n
                expn = collect(multiexponents(length(p.α[i]),k))
                for exponent in expn
                    u_i = 1
                    for (j,e) in enumerate(exponent)
                        u_i*=p.α[i][j]^e
                    end
                    if k > 0
                        push!(α,u_i)
                    end
                    println(u_i)
                    push!(u,p.γ[i]*u_i)
                end
            end
        end
    end 

    println("u:")
    display(u)

    C = Array{Sym,2}(undef,(length(Cu),length(u)))
    
    println(collect(p.γ[i]=>0 for i in 1:length(p.γ)))
    println(collect(α[i]=>0 for i in 1:length(α)))
  
    println(Cu)
    for i in 1:length(Cu)
       
        C[i,1] = Cu[i]((p.γ[k]=>0 for k in 1:length(p.γ))...)
        for j in 2:length(u)
           
            C[i,j] = Cu[i].coeff(u[j])
        end
    end

    println("C:")
    println(C)
    Ca =Array{Sym,2}( C.subs([α[i]=>0 for i in 1:length(α)]))
    println("Ca:")
    println(Ca)

    return Ca,u
end

#Compute constraint from error expression 
function getConstraints(err,n,d,p)

    #Get error of n-th order
    ϵₙ = 0
    if n > 0
        ϵₙ = err.coeff(δ^n)
    else
        ϵₙ = err(δ=>0)
    end

    println("$n-th order error contribution:")
    display(ϵₙ)

    #Factorize the partial derivatives of f
    multiindex = collect(multiexponents(d,n))
    Cu = []
    for idx in multiindex
        df = partial_derivative(idx)
        println(df)
        Cu_i = ϵₙ.coeff(df)
        println(Cu_i)
        if Cu_i != 0
            push!(Cu,Cu_i)
        end
    end

    #display(Cu)

    C,u = getCoeffMatrix(Cu,n,p)

    C_r = rref(C)
    println("Cᵣ:")
    display(C_r)

    C_r = C_r[vec([any(C_r[i,:] .!= 0) for i in 1:size(C_r,1)]),:]

    return C, C_r,u 
end

function getOptimizationEquations(d,s,pQ)
    numParameters = sum(pQ.p.n)

    order=0;
    numEquations = 0
    constraints = []
    objective = 1.0
    high_order = 0

    while numEquations<=numParameters && order<20

        #f = nth_multitaylor(order,d)
        f = nth_deegree_poly(order,d)
        err = error(d,pQ.x,pQ.w,f,s)
        C,C_r,u = getConstraints(err,order,d,pQ.p)
        
        if length(C_r) > 0
            println("Constraints $order-th order:")
            display(C_r*u)
            loc_c = C_r*u

            numEquations+=length(loc_c)

            if numEquations>numParameters
                objective = (C,u) #expand(sum((C*u).^2))
            else
                push!(constraints,loc_c)
                high_order = order
            end

        end

        println("#Eq./#Para.: ",numEquations,"/",numParameters)
        
        order+=1
    end

    return objective, constraints, high_order
end


