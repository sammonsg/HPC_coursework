import matplotlib.pyplot as plt
import numpy as np
import matplotlib

# Constants
E = 210e9
w = 0.100
h = 0.120
I = (w * h ** 3) / 12
Fy = 1000
F_conc = 1000
L = 10

# Calculate theoretical deflections
x = np.linspace(0, L, 101)
u = np.zeros(101)
u = (F_conc * x * x * (L - x) ** 2) / (24 * E * I) + (Fy * x * x * (3 * L - 4 * x)) / (48 * E * I)


# Mirror solution
for x in xrange(1, 50):
    # u3[51 - x] = u[51 - x]
    u[x + 51] = u[50 - x]


# Concentrated Load
plt.title('Deflection Task 1')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.ylabel('Deflection (m)')
plt.xlabel('Span Position (m)')
plt.grid(b=True, which='major', color='k', linestyle='--')
crs = open('./output/data/u_solution.txt', 'r')
D = np.loadtxt(crs)
x2 = np.linspace(0, 10, len(D))
pos = np.linspace(0, 10, 101)
plt.plot(x2, D[:,1], 'k', lw=2, label='Computed')
plt.plot(pos, u, 'g--', lw=1.25, label='Analytical')
plt.legend(bbox_to_anchor=(.5, 0), loc=8, borderaxespad=0., title="Solution")
plt.savefig('.task1', bbox_inches='tight')
