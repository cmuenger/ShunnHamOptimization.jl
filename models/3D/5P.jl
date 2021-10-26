function objective(α1, γ1, γ2)
    1.0
end

function objective_grad1(α1, γ1, γ2)
    0
end

function objective_grad2(α1, γ1, γ2)
    0
end

function objective_grad3(α1, γ1, γ2)
    0
end

objective_grad = [objective_grad1,objective_grad2,objective_grad3]

function constraint1(α1, γ1, γ2)
    (-γ1 - 4γ2) + 1
end

function constraint1_grad1(α1, γ1, γ2)
    0
end

function constraint1_grad2(α1, γ1, γ2)
    -1
end

function constraint1_grad3(α1, γ1, γ2)
    -4
end

constraint1_grad = [constraint1_grad1,constraint1_grad2,constraint1_grad3]

function constraint2(α1, γ1, γ2)
    (((-320 * α1 .^ 2) .* γ2 + (160α1) .* γ2) - 20γ2) + 1
end

function constraint2_grad1(α1, γ1, γ2)
    (-640α1) .* γ2 + 160γ2
end

function constraint2_grad2(α1, γ1, γ2)
    0
end

function constraint2_grad3(α1, γ1, γ2)
    (-320 * α1 .^ 2 + 160α1) - 20
end

constraint2_grad = [constraint2_grad1,constraint2_grad2,constraint2_grad3]

function constraint3(α1, γ1, γ2)
    ((((3840 * α1 .^ 3) .* γ2 - (2880 * α1 .^ 2) .* γ2) + (720α1) .* γ2) - 60γ2) + 1
end

function constraint3_grad1(α1, γ1, γ2)
    ((11520 * α1 .^ 2) .* γ2 - (5760α1) .* γ2) + 720γ2
end

function constraint3_grad2(α1, γ1, γ2)
    0
end

function constraint3_grad3(α1, γ1, γ2)
    ((3840 * α1 .^ 3 - 2880 * α1 .^ 2) + 720α1) - 60
end

constraint3_grad = [constraint3_grad1,constraint3_grad2,constraint3_grad3]

P5_model_j = JuliaModel(FunctionGradient(objective,objective_grad),
	[FunctionGradient(constraint1,constraint1_grad),
	 FunctionGradient(constraint2,constraint2_grad),
	 FunctionGradient(constraint3,constraint3_grad)],
	[0.3, 0.2, 0.2],
	Parameter((1,2))
	)
 
