#!/usr/bin/env julia
# using Revise; using RHEOS; using NLopt; using BenchmarkTools

time_sim = timeline(t_start = 0, t_end = 100, step = 0.1)
load_sim = strainfunction(time_sim, ramp(offset = 0.0, gradient = 1.0));
param_sim = (η = 2.0, kᵦ=0.5, kᵧ=1.0)
model_sim = RheoModel(SLS_Zener, param_sim)

data = modelpredict(load_sim, model_sim)

deriv = RHEOS.derivBD

σ = data.σ
σdot = deriv(σ, data.t)
ϵ = data.ϵ
ϵdot = deriv(ϵ, data.t)

function differentialCost(σ, σdot, ϵ, ϵdot, params)

    η, kᵦ, kᵧ = params

    cost =  sum(((η + η*kᵧ/kᵦ)*ϵdot + kᵧ*ϵ - σ - (η/kᵦ)*σdot).^2)

end

# fitting test to show it actually works
# p0 = [400.0, 0.25, 0.0007]
# opt = Opt(:LN_SBPLX, length(p0))
# xtol_rel!(opt, 1e-10)
# min_objective!(opt, (params, grad) -> differentialCost(σ, σdot, ϵ, ϵdot, params))
# (minf, minx, ret) = NLopt.optimize(opt, p0)
# println(minx)

# benchmark test to compare speeds
dt = data.t[2] - data.t[1];
vals = values(param_sim)

println("Boltzmann cost function:")
@btime RHEOS.obj_const_nonsing($vals, 0, $model_sim._Ga, $data.t, $dt, $ϵdot, $σ) # 110.596 μs

println("\nDifferential Equation cost function:")
@btime differentialCost($σ, $σdot, $ϵ, $ϵdot, $vals) # 5.304 μs

println("")