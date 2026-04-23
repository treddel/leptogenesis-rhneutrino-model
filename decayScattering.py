import numpy as np
from scipy.special import kn
import matplotlib.pyplot as plt


K = 10**1
K_s = 1
a5 = (8* np.pi**2)/(9*np.log(10**5))
a1 = (8* np.pi**2)/(9*np.log(10))
z_span = [10**-4, 10**2]
z = np.linspace(z_span[0], z_span[1], 160764)
D = K*z*(kn(1,z)/kn(2,z))
DS5 = K_s * (1 + np.log(10**5) * z**2 * np.log(1 + a5/z))
DS1 = K_s * (1 + np.log(10) * z**2 * np.log(1 + a1/z))

plt.plot(z, D, 'b-', label='D', linewidth=0.5)
plt.plot(z, DS5, 'r--', label='D+S $( \dfrac{M_h}{M_1} = 10^{-5})$', linewidth=0.5)  # Mh/M1 = 10^-5
plt.plot(z, DS1, 'r-.', label='D+S $( \dfrac{M_h}{M_1} = 10^{-1})$', linewidth=0.5)  # Mh/M1 = 10^-1
plt.xlabel('z = M$_1$ / T')
plt.ylabel('')
plt.xscale('log')
plt.yscale('log')
plt.xlim([10**-2, 10**2])
plt.ylim([10**-2, 10**3])
plt.show()

