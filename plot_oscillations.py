import matplotlib.pyplot as plt
import numpy as np
import matplotlib


# Concentrated Load
plt.title('Deflection Task 2 dt = 0.00005s')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.ylabel('Deflection (m)')
plt.xlabel('Time (s)')
# plt.grid(b=True, which='major', color='k', linestyle='--')
crs = open('./output/data/u_centre_log.txt', 'r')
D = np.loadtxt(crs)
t = np.linspace(0, 5, len(D))
pos = np.linspace(0, 10, 101)
plt.plot(t, D, 'k', lw=2, label='Computed')
# plt.plot(pos, u, 'r--', lw=1.25, label='Analytical')
plt.savefig('.task2', bbox_inches='tight')
