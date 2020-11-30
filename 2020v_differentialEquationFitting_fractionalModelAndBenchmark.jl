#!/usr/bin/env julia
# using Revise; using RHEOS; using NLopt; using DSP; using SpecialFunctions; using BenchmarkTools

# fractSLS_predicted = modelpredict(data, FractionalSLS([2.0, 0.5, 0.5, 0.7]), :G)
# cₐ, a, kᵦ, kᵧ = params

deriv = RHEOS.derivBD

data = loaddata("fractsls_predicted.jld2")
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

function differentialCost(σ, σderiv, σdoublederiv, ϵ, ϵderiv, ϵdoublederiv, t, dt, params, grad)

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
# min_objective!(opt, (params, grad) -> differentialCost(σ, σderiv, σdoublederiv, ϵ, ϵderiv, ϵdoublederiv, data.t, dt, params, grad))
# (minf, minx, ret) = NLopt.optimize(opt, p0)
# println(minx)

# BENCHMARKS
@btime differentialCost($σ, $σderiv, $σdoublederiv, $ϵ, $ϵderiv, $ϵdoublederiv, $data.t, $dt, $p0, 1.0)

testmodel = FractionalSLS()
testmodulus = testmodel.G
@btime RHEOS.obj_const_nonsing($p0, [1.0], $testmodulus, $t, $dt, $ϵderiv, $σ)


