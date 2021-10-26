#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include <nlopt.hpp>

#include "model.hpp"

#include "models/2D/P28.hpp"

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) 
{
    s.put('[');
    char comma[3] = {'\0', ' ', '\0'};
    for (const auto& e : v) {
        s << comma << e;
        comma[0] = ',';
    }
    return s << ']';
}



double myfunction(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    FunctionGradient *d = reinterpret_cast<FunctionGradient*>(data);
    if (!grad.empty()) {
        for(int i=0; i<grad.size(); i++) {
            grad[i] = d->grad[i](x);
        }
    }
    return d->f(x);
}

void test_constraints(Model& model, const std::vector<double>& x) {
    std::cout<<"Constraints: ";
    for(const auto& c : model.constraints)
    {
        std::cout<<c.f(x)<<", ";
    }
    std::cout<<std::endl;
}

void run_opt_sequential(Model& model) {
    int n1 = model.p.first;
    int n2 = model.p.second;

    std::vector<double> lb(n1+n2);
    std::vector<double> ub(n1+n2);
    for(int i = 0; i<n1; i++) {
        lb[i] = 0.0;
        ub[i] = 1.0;
    }
    for(int i = n1; i<n1+n2; i++) {
        lb[i] = 0.0;
        ub[i] = HUGE_VAL;
    }
    //Start with random init
    std::vector<double> x = init_rand(model);
    std::cout<<"Init: "<<x<<std::endl;
    //Optimization constraints first
    for( int i=1; i<model.constraints.size(); i++) {
        std::cout<<i<<" constraints"<<std::endl;
        //Build new objective function
        const auto& objective_f = [model,i] (const std::vector<double>& x) -> double {
            return std::pow(model.constraints[i].f(x),2);
        };
        std::vector<func_t> objective_grad;
        for(int j=0; j<model.objective.grad.size(); j++)
        {
            const auto f = [model,i,j](const std::vector<double>& x) -> double {
                return 2*model.constraints[i].f(x)*model.constraints[i].grad[j](x);
            };

            objective_grad.push_back(f);
        }

        std::vector<FunctionGradient> constraints;
        for(int j=0; j<i; j++)
        {
            constraints.push_back(model.constraints[j]);
        }

        Model model_c{{objective_f,objective_grad},constraints,x,model.p};

    
        //Run optimization
        nlopt::opt opt(nlopt::LD_SLSQP, n1+n2);
    
        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(ub);

        opt.set_min_objective(myfunction, &model_c.objective);

        for(FunctionGradient& c : model_c.constraints) {
            opt.add_equality_constraint(myfunction, &c);
        }
  
        opt.set_xtol_rel(1e-9);
        opt.set_maxeval(1000);
    
        double minf;
        
        try{
            nlopt::result result = opt.optimize(x, minf);
            //std::cout << "found minimum at f("<<x<<") = "<<std::setprecision(10) << minf << std::endl;
        }
        catch(std::exception &e) {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }
    }

    std::cout<<"objective"<<std::endl;
     //Run optimization
    nlopt::opt opt(nlopt::LD_SLSQP, n1+n2);
    
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    opt.set_min_objective(myfunction, &model.objective);

    for(FunctionGradient& c : model.constraints) {
        opt.add_equality_constraint(myfunction, &c);
    }
  
    opt.set_xtol_rel(1e-8);
    opt.set_maxeval(1000);

    double minf;

    try{
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "found minimum at f("<<x<<") = "<<std::setprecision(17) << minf << std::endl;
        test_constraints(model,x);
        std::cout<<std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }


}

void run_opt(Model& model) {
    int n1 = model.p.first;
    int n2 = model.p.second;

    nlopt::opt opt(nlopt::LD_SLSQP, n1+n2);
    std::vector<double> lb(n1+n2);
    std::vector<double> ub(n1+n2);
    for(int i = 0; i<n1; i++) {
        lb[i] = 0.0;
        ub[i] = 1.0;
    }
    for(int i = n1; i<n1+n2; i++) {
        lb[i] = 0.0;
        ub[i] = HUGE_VAL;
    }


    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    opt.set_min_objective(myfunction, &model.objective);

    for(FunctionGradient& c : model.constraints) {
         opt.add_equality_constraint(myfunction, &c);
    }
  
    opt.set_xtol_rel(1e-8);
    opt.set_maxeval(1000);
    
    std::vector<double> x = init_rand(model);
    std::cout<<"Init: "<<x<<std::endl;


    double minf;

    try{
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "found minimum at f(";
        for(double xi : x) {
            std::cout<< xi << ",";
        }
        std::cout<<") = "
            << std::setprecision(10) << minf << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
}


int main()
{
    Model model = P28_model_c;
    

    int M = 50;
    for(int i=0; i<M; i++) {
        run_opt_sequential(model);
    }
    
    
    return 0;
}

