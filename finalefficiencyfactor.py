import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt


y_span = [10**-5, 10**1]
x_span = [10**-3, 10**3]

K = np.linspace(x_span[0], x_span[1], 10000)
zb = 1 + 1/2 * np.log(1 + (np.pi * K**2)/1024 * np.log((3125*np.pi*K**2)/1024)**5)
kappa = 2/(zb*K) * (1-np.exp(-1/2 * zb*K))

plt.plot(K, kappa, 'r-', linewidth=0.5)
#  plt.legend()
plt.xlabel('$K$')
plt.ylabel('$\kappa_f$')
#  plt.title("Washout ")
plt.xscale('log')
plt.yscale('log')
plt.xlim([10**-3, 10**3])
plt.ylim([10**-3, 10**1])
plt.show()
