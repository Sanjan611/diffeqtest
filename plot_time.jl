using PyPlot

# PyPlot.clf()

times = [303.846, 135.059, 75.768, 33.558, 18.952, 12.130, 8.366, 6.125]
costs = [0.000144, 0.000196, 0.000239, 0.000312, 0.000375, 0.000439, 0.000512, 0.000595]
num_points = [1001, 667, 501, 334, 251, 201, 167, 143]

plot(num_points, times)
xlabel("Number of time stamps")
ylabel("Time taken (milliseconds)")
savefig("images/num_times.eps")
display(gcf())

# PyPlot.clf()

plot(num_points, costs)
xlabel("Number of time stamps")
ylabel("Cost")
savefig("images/num_costs.eps")