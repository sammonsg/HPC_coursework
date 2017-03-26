import numpy as np
import matplotlib.pyplot as plt


# Concentrated Load
plt.title('Oscillation on explicit scheme, dt = 0.0005s')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.ylabel('Amplitude of oscillations (m)')
plt.xlabel('Force application time (s)')
plt.grid(b=True, which='major', color='k', linestyle='--')

# Fetch data
crs = open('./output/data/minmax_exp.txt', 'r')
data = np.loadtxt(crs)
plt.loglog(data[:,1], data[:,2], 'k', lw=2, label='Computed')

# Save as...
plt.savefig('task2fff', bbox_inches='tight')
