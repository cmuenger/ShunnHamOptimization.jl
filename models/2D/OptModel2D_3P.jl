P3_s = Parameter((1,1))

objective3_s = (Array{Sym,2}(ones(1,1)),Array{Sym,1}([1.0]))
constraints3_s =Array{Sym}(undef,2)
constraints3_s[1] = -3*γ1 + 1
constraints3_s[2] = -108*α1^2*γ1 + 72*α1*γ1 - 12*γ1 + 1

init3_s = [0.211324865405187, 1/3]

P3_model_s = OptModel(objective3_s, constraints3_s, init3_s, P3_s)