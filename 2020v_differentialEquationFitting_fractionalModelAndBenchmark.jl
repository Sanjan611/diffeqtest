#!/usr/bin/env julia
# using Revise; using RHEOS; using NLopt; using BenchmarkTools; using DSP; using SpecialFunctions

time_sim = timeline(t_start = 0, t_end = 100, step = 0.1)
load_sim = strainfunction(time_sim, ramp(offset = 0.0, gradient = 1.0));
param_sim = (cₐ = 2.0, a = 0.5, kᵦ = 0.5, kᵧ = 0.7)
model_sim = RheoModel(FractSLS_Zener, param_sim)

deriv = RHEOS.derivBD

t = data.t
dt = t[2] - t[1]
σ = data.σ
σderiv = deriv(σ, t)
σdoublederiv = deriv(σderiv, t)
ϵ = data.ϵ
ϵderiv = deriv(ϵ, t)
ϵdoublederiv = deriv(ϵderiv, t)

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
# p0 = [3.0, 0.2, 5.0, 1.0]
# opt = Opt(:LN_SBPLX, length(p0))
# lower_bounds!(opt, [0.0, 0.0, 0.0, 0.0])
# upper_bounds!(opt, [Inf, 1.0, Inf, Inf])
# xtol_rel!(opt, 1e-5)
# min_objective!(opt, (params, grad) -> differentialCost(σ, σderiv, σdoublederiv, ϵ, ϵderiv, ϵdoublederiv, data.t, dt, params))
# (minf, minx, ret) = NLopt.optimize(opt, p0)
# println(minx)

# BENCHMARKS
@btime differentialCost($σ, $σderiv, $σdoublederiv, $ϵ, $ϵderiv, $ϵdoublederiv, $data.t, $dt, $values(param_sim))

testmodulus = model_sim._Ga
dummygrad = [1.0]
@btime RHEOS.obj_const_nonsing($values(param_sim), $dummygrad, $testmodulus, $t, $dt, $ϵderiv, $σ)



