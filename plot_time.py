import matplotlib.pyplot as plt

times = [303.846, 135.059, 75.768, 33.558, 18.952, 12.130, 8.366, 6.125]
costs = [0.144, 0.196, 0.239, 0.312, 0.375, 0.439, 0.512, 0.595]
num_points = [1001, 667, 501, 334, 251, 201, 167, 143]

plt.plot(num_points, times)
plt.scatter(num_points, times)
plt.xlabel("Number of time points")
plt.ylabel("Time taken (milliseconds)")
plt.tight_layout()
plt.grid(True)
plt.savefig("num_times.eps")
plt.show()

plt.plot(num_points, costs)
plt.scatter(num_points, costs)
plt.xlabel("Number of time points")
plt.ylabel("Cost (10^(-3))")
plt.tight_layout()
plt.grid(True)
plt.savefig("num_costs.eps")
plt.show()


#
#  using PyPlot

# # PyPlot.clf()

# times = [303.846, 135.059, 75.768, 33.558, 18.952, 12.130, 8.366, 6.125]
# costs = [0.000144, 0.000196, 0.000239, 0.000312, 0.000375, 0.000439, 0.000512, 0.000595]
# num_points = [1001, 667, 501, 334, 251, 201, 167, 143]

# plot(num_points, times)
# xlabel("Number of time stamps")
# ylabel("Time taken (milliseconds)")
# savefig("images/num_times.eps")
# display(gcf())

# # PyPlot.clf()

# plot(num_points, costs
# xlabel("Number of time stamps")
# ylabel("Cost")
# savefig("images/num_costs.eps")