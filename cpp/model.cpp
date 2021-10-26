#include "model.hpp"

std::vector<double> perturbation(const Model& model) {
    std::random_device rd;  
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform_dis(-1.0, 1.0);

    std::vector<double> dx(model.p.first+model.p.second,0.0);
    for( int i =0; i < model.p.first; i++) {
        dx[i] = uniform_dis(gen);
    }

    return dx;
}

std::vector<double> init_rand(const Model& model) { 
    std::vector<double> x = model.init;
    std::vector<double> dx = perturbation(model);
    for(int i=0; i<x.size(); i++) {
        x[i] += 0.8*x[i]*dx[i];
    }

    return x;
}