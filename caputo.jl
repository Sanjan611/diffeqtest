using SpecialFunctions
using RHEOS
using NLopt; using BenchmarkTools;using DSP

################## L1 operator ##################################

function a(j, α)
    return (j+1)^(1-α) - j^(1-α)
end

function L1(f, t, Δt, α)

    N = size(t,1)
    l1 = zeros(Float64, N)
    g = gamma(2-α)
    
    for k in 1:N-1
        l1_k = f[k+1]
        for j in 0:k-1
            l1_k = l1_k - (a(k-j-1, α) - a(k-j, α))*f[j+1]
        end
        l1_k = (l1_k - f[1]*a(k-1, α))*(Δt^(-α))/g
        # l1_k = l1_k*(dt^(-α))/g
        l1[k] = l1_k
    end

    return l1
end

##################### L1-2 Operator ############################

function b(j, α)
    # j >= 0
    return (1/(2-α))*((j+1)^(2-α) - j^(2-α)) - 0.5*((j+1)^(1-α) + j^(1-α))
end

function δt(f, k, Δt)
    # the paper implementation used δₜf_{k-1/2}  
    # modified a bit for code purposes
    return (f[k]-f[k-1])/Δt
end

function δ2t(f, k, Δt)
    return (δt(f,k+1,Δt)-δt(f,k,Δt))/Δt
end

function L12(f, t, Δt, α)

    N = size(t,1)
    l12 = zeros(Float64, N)
    L1_ = L1(f, t, Δt, α)
    g = gamma(2-α)

    for k in 1:N
        l12_k = 0
        if k>2 # only then k-j will be valid
            for j in 2:k-1
                l12_k = l12_k + b(k-j, α)*δ2t(f, (j-1)+1, Δt) # extra +1 to handle Julia indexing vs paper indexing
            end
        end
        l12_k = l12_k*(Δt^(2-α))/g
        l12[k] = l12_k + L1_[k]

    end

    return l12
    
end




function fractionalDiff(sigderiv, sigdoublederiv, t, dt, α)

    I1 = sigderiv[1]*t.^(1-α)
    I2 = conv(t.^(1-α), sigdoublederiv)[1:length(t)]*dt
    return (I1 + I2)/gamma(2 - α)

end
