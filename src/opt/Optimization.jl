function build_function_julia(f)

    f_o = function (x::Vector, grad::Vector)  
                    if length(grad) > 0
                        for i in 1:length(grad)
                            grad[i] =f.grad[i](x...)
                        end
                    end
                    return f.f(x...)
                end
    return f_o
end


function build_objective_2(f, ϕ)

    #Cu = f[1]*f[2]
    #objective =sum(Cu.^2)
    fC = lambdify(f[1],[])
    fu = lambdify(f[2],ϕ)
    #fn = lambdify(objective,ϕ)

    gradu = []
    for i in 1:length(ϕ)
        push!(gradu,lambdify(diff.(f[2],ϕ[i]),ϕ))
        # push!(gradu,lambdify(diff(fn,ϕ[i]),ϕ))
    end

    f_o = function (x::Vector, grad::Vector)  
                    Cu = fC()*fu(x...)
                    if length(grad) > 0
                        for i in 1:length(grad)
                            Cgu = fC()*gradu[i](x...)
                            grad[i] = 2*dot(Cu,Cgu)
                        end
                    end
                    return dot(Cu,Cu)
                end
    return f_o
end

function build_objective_3(f, ϕ)

    Cu = f[1]*f[2]
    objective = sum(Cu.^2)
   
    fn = lambdify(objective,ϕ)

    gradu = []
    for i in 1:length(ϕ)
        push!(gradu,lambdify(diff(objective,ϕ[i]),ϕ))
    end

    f_o = function (x::Vector, grad::Vector)  
                    if length(grad) > 0
                        for i in 1:length(grad)
                            grad[i] = gradu[i](x...)
                        end
                    end
                    return fn(x...)
                end
    return f_o
end

function build_function_2(f, ϕ)

    fn = lambdify(f,ϕ)

    gradn = []
    for i in 1:length(ϕ)
        push!(gradn,lambdify(diff(f,ϕ[i]),ϕ))
    end

    f_c = function (x::Vector, grad::Vector)  
                    if length(grad) > 0
                        for i in 1:length(grad)
                            grad[i] = gradn[i](x...)
                        end
                    end
                    return fn(x...)
                end
    return f_c
end


function compute_constraint(f, ϕ,x)
    subphi = (ϕ[j]=>x[j] for j in 1:length(ϕ))
    return f((subphi...))
end

function test_constraints(model::OptModel,x)

    ϕ = get_phi(model)
    r = []
    for c in model.constraints
        push!(r,compute_constraint(c, ϕ, x))  
    end

    return r
end


function test_constraints(model::JuliaModel,x)

    r = []
    for c in model.constraints
        f_c = build_function_julia(c)
        push!(r,f_c(x,[]))  
    end

    return r
end

function find_init(model, m =250)

    best_r = Inf
    best_x = []
    for i in 1:m
         dx=perturbation(model)
            
        x = model.init +1*model.init.*dx
    
        r = test_constraints(model,x)
        if norm(r) < best_r
            best_r = norm(r)
            best_x = x
        end
    end
    return (best_x, best_r)
end


function id_init(model, m =250)
    return (model.init, 0.0)
end



function run_opt_sequential(model::JuliaModel) 
    
    (x,r) = find_init(model,1)
    new_init = x
    println(new_init)
    for i in 2:length(model.constraints)
        println("$(i-1) constraints")
        #new objective
        new_objective_f = (x...) -> model.constraints[i].f(x...)^2
        new_objective_grad = []
        for g in model.constraints[i].grad
            new_grad_i = (x...) -> 2*model.constraints[i].f(x...)*g(x...)
            push!(new_objective_grad,new_grad_i)
        end

        #new constraint
        new_constraint = model.constraints[1:i-1]

        c_model = JuliaModel(FunctionGradient(new_objective_f,new_objective_grad),new_constraint,new_init,model.p)
      

        (g,x) = run_opt(c_model,1,1,id_init)

        new_init = x
    end

    o_model = JuliaModel(model.objective,model.constraints,new_init,model.p)
    println("Objective")
    return run_opt(o_model,1,1,id_init)
    
end

function run_opt(model::OptModel,n=1, m=1, init=find_init,  maxeval = 500)

    ϕ = get_phi(model)

    opt = Opt(:LD_SLSQP, numParameters(model))
    opt.lower_bounds = lower_bounds(model)
    opt.upper_bounds = upper_bounds(model)
    opt.xtol_rel = 1e-8
    opt.maxeval = maxeval
 

    opt.min_objective =  build_objective_2(model.objective,ϕ)


    for c in model.constraints
        f_c = build_function_2(c,ϕ)
        equality_constraint!(opt, (x,g) -> f_c(x,g))
    end

    minf_g = Inf;
    minx_g = []

    for i in 1:n
        (x,r) = init(model,m)
        println("init $x ($r)")
        (minf,minx,ret) = optimize(opt, x)
        numevals = opt.numevals # the number of function evaluations
        println("got $minf at $minx after $numevals iterations (returned $ret)")
        #println("sum wts: ", 3*minx[8]+3*minx[9]+3*minx[10]+6*minx[11]+6*minx[12])
        r = test_constraints(model, minx)
        println("constaints error $r")

        if minf < minf_g
            minf_g = minf
            minx_g = minx
        end
    end

    return (minf_g,minx_g)
end


function run_opt(model::JuliaModel,  n=1, m=1, init=find_init, maxeval=10000)

  
    opt = Opt(:LD_SLSQP, numParameters(model))
    opt.lower_bounds = lower_bounds(model)
    opt.upper_bounds = upper_bounds(model)
    opt.xtol_rel = 1e-8
    #opt.ftol_rel = 1e-9
    opt.maxeval = maxeval
 

    opt.min_objective =  build_function_julia(model.objective)


    for c in model.constraints
        f_c =  build_function_julia(c)
        equality_constraint!(opt, (x,g) -> f_c(x,g))
    end

    minf_g = Inf;
    minx_g = []

    for i in 1:n
        (x,r) = init(model,m)
        #println("init $x ($r)")
        (minf,minx,ret) = optimize(opt, x)
        numevals = opt.numevals # the number of function evaluations
        #println("got $minf at $minx after $numevals iterations (returned $ret)")
        #println("sum wts: ", 3*minx[8]+3*minx[9]+3*minx[10]+6*minx[11]+6*minx[12])
        r = test_constraints(model, minx)
        #println("constaints error $r")

        if minf < minf_g
            minf_g = minf
            minx_g = minx
        end
    end

    return (minf_g,minx_g)
end

function find_opt(model)
    minf_g = Inf;
    minx_g = []

    for i in 1:5000
    
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
