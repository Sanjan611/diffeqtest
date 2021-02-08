include("caputo.jl")

using DataFrames
using CSV

df = DataFrame(α = Float32[], E1 = Float32[], E1r = Float32[], E2 = Float32[], E2r = Float32[])

function test(t, α = 0.5)
    f = t.^(4+α)
    return f
end

function test_ans(α)
    return gamma(5+α)/24
end


αs = [0.9, 0.5, 0.1]
Δts = [1/10, 1/20, 1/40, 1/80, 1/160, 1/320, 1/640, 1/1280, 1/2560, 1/5120]

for α in αs
    println(α)
    for (i, Δt) in enumerate(Δts)
        println("--- ", Δt)
        t = collect(0:Δt:2)
        f = test(t, α)
        f1 = L1(f, t, Δt, α)
        f2 = L12(f, t, Δt, α)
        error1 = abs(test_ans(α) - f1[findfirst(isequal(1), t)])
        error2 = abs(test_ans(α) - f2[findfirst(isequal(1), t)])
        # println(i, " -- ", error1, " -- ", error2)
        push!(df, (α, error1, 0, error2, 0))
    end
end

for i in 1:size(df,1)-1
    df.E1r[i] = log(2, df.E1[i]/df.E1[i+1])
    df.E2r[i] = log(2, df.E2[i]/df.E2[i+1])
end

CSV.write("test_example.csv", df)
