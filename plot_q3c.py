import numpy as np
import matplotlib.pyplot as plt



plt.title('Oscillation on implicit scheme, dt = 0.005s')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.ylabel('Amplitude of oscillations (m)')
plt.xlabel('Force application time (s)')
plt.grid(b=True, which='major', color='k', linestyle='--')

# Fetch data
crs = open('./output/data/minmax_imp.txt', 'r')
data = np.loadtxt(crs)
plt.loglog(data[:,1], data[:,2], 'k', lw=2, label='Computed')

# Save as...
plt.savefig('task3', bbox_inches='tight')
