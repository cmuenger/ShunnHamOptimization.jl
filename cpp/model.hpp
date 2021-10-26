
#ifndef MODEL_HPP
#define MODEL_HPP

#include <cmath>
#include <functional>
#include <random>
#include <vector>

typedef std::function<double(const std::vector<double>&)> func_t;
struct FunctionGradient {
	 func_t f;
	 std::vector<func_t> grad;
};

struct Model {
	 FunctionGradient objective;
	 std::vector<FunctionGradient> constraints;
	 std::vector<double> init;
	 std::pair<int,int> p;
};

std::vector<double> perturbation(const Model& model); 

std::vector<double> init_rand(const Model& model);

#endif /* MODEL_HPP */