#!/usr/bin/env julia
# using RHEOS
# using NLopt

# sls_predicted = modelpredict(data, SLS([2.0, 0.5, 1.0]), :G)
# η, kᵦ, kᵧ

deriv = RHEOS.derivBD

data = loaddata("sls_predicted.jld2")
σ = data.σ
σdot = deriv(σ, data.t)
ϵ = data.ϵ
ϵdot = deriv(ϵ, data.t)

function differentialCost(σ, σdot, ϵ, ϵdot, params, grad)

    η, kᵦ, kᵧ = params

    cost =  sum(((η + η*kᵧ/kᵦ)*ϵdot + kᵧ*ϵ - σ - (η/kᵦ)*σdot).^2)

end

# original parameters = [2.0, 0.5, 1.0]
p0 = [400.0, 0.25, 0.0007]
opt = Opt(:LN_SBPLX, length(p0))
# lower_bounds!(opt, low_bounds)
# upper_bounds!(opt, hi_bounds)
xtol_rel!(opt, 1e-10)
min_objective!(opt, (params, grad) -> differentialCost(σ, σdot, ϵ, ϵdot, params, grad))
(minf, minx, ret) = NLopt.optimize(opt, p0)
println(minx)