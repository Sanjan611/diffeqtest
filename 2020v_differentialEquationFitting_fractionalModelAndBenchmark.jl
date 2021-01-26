#!/usr/bin/env julia
# using Revise;
using RHEOS; using NLopt; using BenchmarkTools; using DSP; using SpecialFunctions

time_sim = timeline(t_start = 0, t_end = 100, step = 0.1)
load_sim = strainfunction(time_sim, ramp(offset = 0.0, gradient = 1.0));
param_sim = (cₐ = 2.0, a = 0.5, kᵦ = 0.5, kᵧ = 0.7)
model_sim = RheoModel(FractSLS_Zener, param_sim)

data = modelpredict(load_sim, model_sim)

deriv = RHEOS.derivBD

t = data.t
dt = t[2] - t[1]
σ = data.σ
σderiv = deriv(σ, t)
σdoublederiv = RHEOS.doublederivCD(σderiv, t)
ϵ = data.ϵ
ϵderiv = deriv(ϵ, t)
ϵdoublederiv = RHEOS.doublederivCD(ϵderiv, t)

function fractionalDiff(sigderiv, sigdoublederiv, t, dt, α)

    I1 = sigderiv[1]*t.^(1-α)
    I2 = conv(t.^(1-α), sigdoublederiv)[1:length(t)]*dt
    return (I1 + I2)/gamma(2 - α)

end

function differentialCost(σ, σderiv, σdoublederiv, ϵ, ϵderiv, ϵdoublederiv, t, dt, params)

    cₐ, a, kᵦ, kᵧ = params

    ϵfdot = fractionalDiff(ϵderiv, ϵdoublederiv, t, dt, a)
    σfdot = fractionalDiff(σderiv, σdoublederiv, t, dt, a)

    # print(params, "\n")
    cost =  sum(((cₐ + cₐ*kᵧ/kᵦ).*ϵfdot .+ kᵧ*ϵ .- σ .- (cₐ/kᵦ).*σfdot).^2)
 
end

# original parameters = [2.0, 0.5, 0.5, 0.7]
# fitting result is drastically afftected by double derivative method used,
# makes sense, particularly as ramp load doesn't have much instantaneous 'information'
# as well.
# p0 = [3.0, 0.6, 0.6, 1.0]
# opt = Opt(:LN_SBPLX, length(p0))
# lower_bounds!(opt, [0.0, 0.0, 0.0, 0.0])
# upper_bounds!(opt, [100, 1.0, 100, 100])
# xtol_rel!(opt, 1e-5)
# min_objective!(opt, (params, grad) -> differentialCost(σ, σderiv, σdoublederiv, ϵ, ϵderiv, ϵdoublederiv, data.t, dt, params))
# (minf, minx, ret) = NLopt.optimize(opt, p0)
# println(minx)

# BENCHMARKS
vals = [i for i in values(param_sim)]

testmodulus = RHEOS._Ga(FractSLS_Zener)
dummygrad = [1.0]

# println("Boltzmann cost function:")
@btime RHEOS.obj_const_nonsing($vals, $dummygrad, $testmodulus, $t, $dt, $ϵderiv, $σ) # 15.858 ms

println("\nDifferential Equation cost function:")
# differentialCost(σ, σderiv, σdoublederiv, ϵ, ϵderiv, ϵdoublederiv, data.t, dt, vals)
@btime differentialCost($σ, $σderiv, $σdoublederiv, $ϵ, $ϵderiv, $ϵdoublederiv, $data.t, $dt, $vals) # 240.779 μs

println("")
