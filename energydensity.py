import numpy as np
import matplotlib.pyplot as plt

T_span = [10**-2, 10**2]
T = np.linspace(T_span[0], T_span[1], 10000)
g = 106.75
rho_r = ((np.pi)**2)/30 * g * T**-4
rho_phi = ((np.pi)**2)/30 * g * T**-(4+2)

plt.plot(T, rho_r, 'k-', label='$ρ_r$', linewidth=0.5)
plt.plot(T, rho_phi, 'r-', label='$ρ_φ$', linewidth=0.5)
plt.legend()
plt.xlabel('$z=M_1 / T$')
plt.ylabel('$ρ(z)$')
plt.xscale('log')
plt.yscale('log')
plt.xlim([10**-2, 10**2])
plt.ylim([10**-6, 10**12])
plt.show()
