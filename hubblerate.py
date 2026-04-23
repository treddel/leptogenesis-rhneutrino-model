import numpy as np
import matplotlib.pyplot as plt

z_span = [10**-4, 10**4]
z = np.linspace(z_span[0], z_span[1], 10000)
g = 106.75
delta = 2
gamma = 50
rho_r = ((np.pi)**2)/30 * g * z**-4
rho_phi = ((np.pi)**2)/30 * g * z**-(4+delta)

H_r = np.sqrt(8*np.pi/3)*np.sqrt(rho_r)
f = np.sqrt(1 + (gamma/z)**delta)
H_mod = H_r*f

plt.plot(z, H_r, 'k-', label='$H_r(z)$', linewidth=0.5)
plt.plot(z, H_mod, 'r-', label='$H^\prime (z)$', linewidth=0.5)
plt.legend()
plt.xlabel('$z=M_1 / T$')
plt.ylabel('$H(z)$')
plt.xscale('log')
plt.yscale('log')
plt.xlim([z_span[0], z_span[1]])
plt.ylim([10**-6, 10**12])
plt.show()
