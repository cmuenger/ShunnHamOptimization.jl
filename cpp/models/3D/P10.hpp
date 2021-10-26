#ifndef  P10_HPP
#define P10_HPP

#include "model.hpp"

double objective(const std::vector<double>& x);

double objective_grad1(const std::vector<double>& x);

double objective_grad2(const std::vector<double>& x);

double objective_grad3(const std::vector<double>& x);

double objective_grad4(const std::vector<double>& x);

double objective_grad5(const std::vector<double>& x);

double objective_grad6(const std::vector<double>& x);

extern std::vector<func_t> objective_grad;

double constraint1(const std::vector<double>& x);

double constraint1_grad1(const std::vector<double>& x);

double constraint1_grad2(const std::vector<double>& x);

double constraint1_grad3(const std::vector<double>& x);

double constraint1_grad4(const std::vector<double>& x);

double constraint1_grad5(const std::vector<double>& x);

double constraint1_grad6(const std::vector<double>& x);

extern std::vector<func_t> constraint1_grad;

double constraint2(const std::vector<double>& x);

double constraint2_grad1(const std::vector<double>& x);

double constraint2_grad2(const std::vector<double>& x);

double constraint2_grad3(const std::vector<double>& x);

double constraint2_grad4(const std::vector<double>& x);

double constraint2_grad5(const std::vector<double>& x);

double constraint2_grad6(const std::vector<double>& x);

extern std::vector<func_t> constraint2_grad;

double constraint3(const std::vector<double>& x);

double constraint3_grad1(const std::vector<double>& x);

double constraint3_grad2(const std::vector<double>& x);

double constraint3_grad3(const std::vector<double>& x);

double constraint3_grad4(const std::vector<double>& x);

double constraint3_grad5(const std::vector<double>& x);

double constraint3_grad6(const std::vector<double>& x);

extern std::vector<func_t> constraint3_grad;

double constraint4(const std::vector<double>& x);

double constraint4_grad1(const std::vector<double>& x);

double constraint4_grad2(const std::vector<double>& x);

double constraint4_grad3(const std::vector<double>& x);

double constraint4_grad4(const std::vector<double>& x);

double constraint4_grad5(const std::vector<double>& x);

double constraint4_grad6(const std::vector<double>& x);

extern std::vector<func_t> constraint4_grad;

double constraint5(const std::vector<double>& x);

double constraint5_grad1(const std::vector<double>& x);

double constraint5_grad2(const std::vector<double>& x);

double constraint5_grad3(const std::vector<double>& x);

double constraint5_grad4(const std::vector<double>& x);

double constraint5_grad5(const std::vector<double>& x);

double constraint5_grad6(const std::vector<double>& x);

extern std::vector<func_t> constraint5_grad;

extern Model P10_model_c;

#endif /* P10_HPP */