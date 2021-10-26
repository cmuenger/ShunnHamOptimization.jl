function objective(α1, α2, α3, γ1, γ2, γ3)
    1.0
end

function objective_grad1(α1, α2, α3, γ1, γ2, γ3)
    0
end

function objective_grad2(α1, α2, α3, γ1, γ2, γ3)
    0
end

function objective_grad3(α1, α2, α3, γ1, γ2, γ3)
    0
end

function objective_grad4(α1, α2, α3, γ1, γ2, γ3)
    0
end

function objective_grad5(α1, α2, α3, γ1, γ2, γ3)
    0
end

function objective_grad6(α1, α2, α3, γ1, γ2, γ3)
    0
end

objective_grad = [objective_grad1,objective_grad2,objective_grad3,objective_grad4,objective_grad5,objective_grad6]

function constraint1(α1, α2, α3, γ1, γ2, γ3)
    ((-4γ1 - 4γ2) - 6γ3) + 1
end

function constraint1_grad1(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint1_grad2(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint1_grad3(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint1_grad4(α1, α2, α3, γ1, γ2, γ3)
    -4
end

function constraint1_grad5(α1, α2, α3, γ1, γ2, γ3)
    -4
end

function constraint1_grad6(α1, α2, α3, γ1, γ2, γ3)
    -6
end

constraint1_grad = [constraint1_grad1,constraint1_grad2,constraint1_grad3,constraint1_grad4,constraint1_grad5,constraint1_grad6]

function constraint2(α1, α2, α3, γ1, γ2, γ3)
    (((((((((-320 * α1 .^ 2) .* γ1 + (160α1) .* γ1) - (320 * α2 .^ 2) .* γ2) + (160α2) .* γ2) - (160 * α3 .^ 2) .* γ3) + (80α3) .* γ3) - 20γ1) - 20γ2) - 10γ3) + 1
end

function constraint2_grad1(α1, α2, α3, γ1, γ2, γ3)
    (-640α1) .* γ1 + 160γ1
end

function constraint2_grad2(α1, α2, α3, γ1, γ2, γ3)
    (-640α2) .* γ2 + 160γ2
end

function constraint2_grad3(α1, α2, α3, γ1, γ2, γ3)
    (-320α3) .* γ3 + 80γ3
end

function constraint2_grad4(α1, α2, α3, γ1, γ2, γ3)
    (-320 * α1 .^ 2 + 160α1) - 20
end

function constraint2_grad5(α1, α2, α3, γ1, γ2, γ3)
    (-320 * α2 .^ 2 + 160α2) - 20
end

function constraint2_grad6(α1, α2, α3, γ1, γ2, γ3)
    (-160 * α3 .^ 2 + 80α3) - 10
end

constraint2_grad = [constraint2_grad1,constraint2_grad2,constraint2_grad3,constraint2_grad4,constraint2_grad5,constraint2_grad6]

function constraint3(α1, α2, α3, γ1, γ2, γ3)
    (((((((3840 * α1 .^ 3) .* γ1 - (2880 * α1 .^ 2) .* γ1) + (720α1) .* γ1 + (3840 * α2 .^ 3) .* γ2) - (2880 * α2 .^ 2) .* γ2) + (720α2) .* γ2) - 60γ1) - 60γ2) + 1
end

function constraint3_grad1(α1, α2, α3, γ1, γ2, γ3)
    ((11520 * α1 .^ 2) .* γ1 - (5760α1) .* γ1) + 720γ1
end

function constraint3_grad2(α1, α2, α3, γ1, γ2, γ3)
    ((11520 * α2 .^ 2) .* γ2 - (5760α2) .* γ2) + 720γ2
end

function constraint3_grad3(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint3_grad4(α1, α2, α3, γ1, γ2, γ3)
    ((3840 * α1 .^ 3 - 2880 * α1 .^ 2) + 720α1) - 60
end

function constraint3_grad5(α1, α2, α3, γ1, γ2, γ3)
    ((3840 * α2 .^ 3 - 2880 * α2 .^ 2) + 720α2) - 60
end

function constraint3_grad6(α1, α2, α3, γ1, γ2, γ3)
    0
end

constraint3_grad = [constraint3_grad1,constraint3_grad2,constraint3_grad3,constraint3_grad4,constraint3_grad5,constraint3_grad6]

function constraint4(α1, α2, α3, γ1, γ2, γ3)
    (((((-13440 * α3 .^ 4) .* γ3 + (13440 * α3 .^ 3) .* γ3) - (5040 * α3 .^ 2) .* γ3) + (840α3) .* γ3) - (105γ3) / 2) + 1
end

function constraint4_grad1(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint4_grad2(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint4_grad3(α1, α2, α3, γ1, γ2, γ3)
    (((-53760 * α3 .^ 3) .* γ3 + (40320 * α3 .^ 2) .* γ3) - (10080α3) .* γ3) + 840γ3
end

function constraint4_grad4(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint4_grad5(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint4_grad6(α1, α2, α3, γ1, γ2, γ3)
    (((-13440 * α3 .^ 4 + 13440 * α3 .^ 3) - 5040 * α3 .^ 2) + 840α3) - 105 / 2
end

constraint4_grad = [constraint4_grad1,constraint4_grad2,constraint4_grad3,constraint4_grad4,constraint4_grad5,constraint4_grad6]

function constraint5(α1, α2, α3, γ1, γ2, γ3)
    ((((((((((((256 * α1 .^ 4) .* γ1 - (256 * α1 .^ 3) .* γ1) + (96 * α1 .^ 2) .* γ1) - (16α1) .* γ1) + (256 * α2 .^ 4) .* γ2) - (256 * α2 .^ 3) .* γ2) + (96 * α2 .^ 2) .* γ2) - (16α2) .* γ2) - (160 * α3 .^ 4) .* γ3) + (160 * α3 .^ 3) .* γ3) - (60 * α3 .^ 2) .* γ3) + (10α3) .* γ3 + γ1 + γ2) - (5γ3) / 8
end

function constraint5_grad1(α1, α2, α3, γ1, γ2, γ3)
    (((1024 * α1 .^ 3) .* γ1 - (768 * α1 .^ 2) .* γ1) + (192α1) .* γ1) - 16γ1
end

function constraint5_grad2(α1, α2, α3, γ1, γ2, γ3)
    (((1024 * α2 .^ 3) .* γ2 - (768 * α2 .^ 2) .* γ2) + (192α2) .* γ2) - 16γ2
end

function constraint5_grad3(α1, α2, α3, γ1, γ2, γ3)
    (((-640 * α3 .^ 3) .* γ3 + (480 * α3 .^ 2) .* γ3) - (120α3) .* γ3) + 10γ3
end

function constraint5_grad4(α1, α2, α3, γ1, γ2, γ3)
    (((256 * α1 .^ 4 - 256 * α1 .^ 3) + 96 * α1 .^ 2) - 16α1) + 1
end

function constraint5_grad5(α1, α2, α3, γ1, γ2, γ3)
    (((256 * α2 .^ 4 - 256 * α2 .^ 3) + 96 * α2 .^ 2) - 16α2) + 1
end

function constraint5_grad6(α1, α2, α3, γ1, γ2, γ3)
    (((-160 * α3 .^ 4 + 160 * α3 .^ 3) - 60 * α3 .^ 2) + 10α3) - 5 / 8
end

constraint5_grad = [constraint5_grad1,constraint5_grad2,constraint5_grad3,constraint5_grad4,constraint5_grad5,constraint5_grad6]

function constraint6(α1, α2, α3, γ1, γ2, γ3)
    (((((((((((143360 * α1 .^ 5) .* γ1 - (179200 * α1 .^ 4) .* γ1) + (89600 * α1 .^ 3) .* γ1) - (22400 * α1 .^ 2) .* γ1) + (2800α1) .* γ1 + (143360 * α2 .^ 5) .* γ2) - (179200 * α2 .^ 4) .* γ2) + (89600 * α2 .^ 3) .* γ2) - (22400 * α2 .^ 2) .* γ2) + (2800α2) .* γ2) - 140γ1) - 140γ2) + 1
end

function constraint6_grad1(α1, α2, α3, γ1, γ2, γ3)
    ((((716800 * α1 .^ 4) .* γ1 - (716800 * α1 .^ 3) .* γ1) + (268800 * α1 .^ 2) .* γ1) - (44800α1) .* γ1) + 2800γ1
end

function constraint6_grad2(α1, α2, α3, γ1, γ2, γ3)
    ((((716800 * α2 .^ 4) .* γ2 - (716800 * α2 .^ 3) .* γ2) + (268800 * α2 .^ 2) .* γ2) - (44800α2) .* γ2) + 2800γ2
end

function constraint6_grad3(α1, α2, α3, γ1, γ2, γ3)
    0
end

function constraint6_grad4(α1, α2, α3, γ1, γ2, γ3)
    ((((143360 * α1 .^ 5 - 179200 * α1 .^ 4) + 89600 * α1 .^ 3) - 22400 * α1 .^ 2) + 2800α1) - 140
end

function constraint6_grad5(α1, α2, α3, γ1, γ2, γ3)
    ((((143360 * α2 .^ 5 - 179200 * α2 .^ 4) + 89600 * α2 .^ 3) - 22400 * α2 .^ 2) + 2800α2) - 140
end

function constraint6_grad6(α1, α2, α3, γ1, γ2, γ3)
    0
end

constraint6_grad = [constraint6_grad1,constraint6_grad2,constraint6_grad3,constraint6_grad4,constraint6_grad5,constraint6_grad6]

P14_model_j = JuliaModel(FunctionGradient(objective,objective_grad),
	[FunctionGradient(constraint1,constraint1_grad),
	 FunctionGradient(constraint2,constraint2_grad),
	 FunctionGradient(constraint3,constraint3_grad),
	 FunctionGradient(constraint4,constraint4_grad),
	 FunctionGradient(constraint5,constraint5_grad),
	 FunctionGradient(constraint6,constraint6_grad)],
	[0.1, 0.2, 0.2, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142],
	Parameter((3,3))
	)
 
