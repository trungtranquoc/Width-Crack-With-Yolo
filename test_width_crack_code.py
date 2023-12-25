from width_crack import *
import time

start_time = time.time()

pol = [[0,5], [4,5], [2,8], [5,6], [5,9], [7,6], [6,3], [7,1], [4,1], [3,0], [1,0], [1,-3], [2,-5], [-2,-5], [-2,-3], [-4,-3], [-2,-1], [-4,1], [-5,-2], [-5,3]]
p = [3,3]

# min_distance, min_x, min_y, vis_pol = optimal_width_crack(p, pol)

# print("Optimal way's result: ", min_distance)

min_distance, min_x, min_y = width_crack(p, pol)

# Print result
print("Result: ", min_distance)

end_time = time.time()
print("Elapse time: ", end_time - start_time)
# plot polygon
pol.append(pol[0])
x, y = zip(*pol)
plt.plot(x, y, color='orange', linewidth=0.5)

# plot p
plt.plot(p[0], p[1], marker='o')

# plot visibility_polygon
# vis_pol.append(vis_pol[0])
# a, b = zip(*vis_pol)
# plt.plot(a, b, color='red', linewidth=0.8)

# plot u,v
plt.plot(min_x[0], min_x[1], marker='o')
plt.plot(min_y[0], min_y[1], marker='o')
u = [min_x[0], min_y[0]]
v = [min_x[1], min_y[1]]
plt.plot(u, v, color='yellow', linewidth=0.5)

plt.show()