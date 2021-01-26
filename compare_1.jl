include("caputo.jl")

using RHEOS
using PyPlot

time_sim = timeline(t_start = 0, t_end = 100, step = 0.15);
# load_sim = strainfunction(time_sim, ramp(offset = 0.0, gradient = 1.0));
# load_sim = strainfunction(time_sim, hstep(offset=5.0, amp=5))
load_sim = strainfunction(time_sim, ramp(offset=0.0, gradient=2.0))  - strainfunction(time_sim, ramp(offset=4, gradient=1.0))
# param_sim = (cₐ = 2.0, a = 0.5, kᵦ = 0.5, kᵧ = 0.7)
param_sim = (cₐ = 2.0, a = 0.6, kᵦ = 0.4, kᵧ = 0.7)
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

# myplot = plot(t, ϵ)
# display(myplot)

α = 0.7 # same as a used in param_sim

function differentialCost1(σ, ϵ, t, dt, params)

    cₐ, a, kᵦ, kᵧ = params

    # calculating fd using Louis' integration by parts approach
    ϵfdot1 = fractionalDiff(ϵderiv, ϵdoublederiv, t, dt, a)
    σfdot1 = fractionalDiff(σderiv, σdoublederiv, t, dt, a)

    cost1 =  sum(((cₐ + cₐ*kᵧ/kᵦ).*ϵfdot1 .+ kᵧ*ϵ .- σ .- (cₐ/kᵦ).*σfdot1).^2)

    return cost1
end

function differentialCost2(σ, ϵ, t, dt, params)

    cₐ, a, kᵦ, kᵧ = params
    α = a

    # calculating fd using L1 discretisation
    ϵfdot2 = L1(ϵ, t, dt, α)
    σfdot2 = L1(σ, t, dt, α)

    cost2 = sum(((cₐ + cₐ*kᵧ/kᵦ).*ϵfdot2 .+ kᵧ*ϵ .- σ .- (cₐ/kᵦ).*σfdot2).^2)/size(t,1)
    return cost2
end

function differentialCost3(σ, ϵ, t, dt, params)
    cₐ, a, kᵦ, kᵧ = params
    α = a

    # calculating fd using L1 discretisation
    ϵfdot3 = L12(ϵ, t, dt, α)
    σfdot3 = L12(σ, t, dt, α)

    cost3 = sum(((cₐ + cₐ*kᵧ/kᵦ).*ϵfdot3 .+ kᵧ*ϵ .- σ .- (cₐ/kᵦ).*σfdot3).^2)/size(t,1)
    return cost3
end

function convIntegralCost(σ, ϵ, t, dt, params)

    cₐ, a, kᵦ, kᵧ = params
    α = a
    data_no_σ = RheoTimeData(ϵ=ϵ, t=t)
    model = RheoModel(FractSLS_Zener, cₐ=cₐ, a=a, kᵦ=kᵦ, kᵧ=kᵧ)
    data_p = modelpredict(data_no_σ, model)

    cost3 = 0
    return cost3
end


# println("Differential Equation cost 1: ", differentialCost1(σ, ϵ, t, dt, param_sim))
# @btime differentialCost1($σ, $ϵ, $t, $dt, $param_sim)

println("Differential Equation cost 2: ", differentialCost2(σ, ϵ, t, dt, param_sim))
@btime differentialCost2($σ, $ϵ, $t, $dt, $param_sim)

println("Differential Equation cost 3: ", differentialCost3(σ, ϵ, t, dt, param_sim))
@btime differentialCost2($σ, $ϵ, $t, $dt, $param_sim)


# as = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
# for a in as
#     param_sim = (cₐ = 2.0, a = a, kᵦ = 0.5, kᵧ = 0.7)
#     println("Differential Equation cost 2: ", differentialCost2(σ, ϵ, t, dt, param_sim))
#     @btime differentialCost2($σ, $ϵ, $t, $dt, $param_sim)
# end

