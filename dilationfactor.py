import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt

z_span = [10**-4, 10**2]
z = np.linspace(z_span[0], z_span[1], 1000)
origin = np.zeros_like(z)
dilation = kn(1, z)/kn(2, z)

plt.plot(z, dilation, 'r', linewidth=0.5)
plt.axhline(y=0.3706, color='k', linewidth=0.5)
plt.axvline(x=1, color='k', linewidth=0.5)
# plt.legend()
plt.xlabel('z = M$_1$ / T')
plt.ylabel('<1/$\gamma$>')
plt.xscale('log')
plt.yscale('log')
plt.xlim([10**-2, 10**2])
plt.ylim([10**-3, 10**(1/2)])

plt.show()

plt.close('all')
