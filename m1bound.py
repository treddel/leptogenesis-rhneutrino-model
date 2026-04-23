import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt

K = 100
gamma = 50
delta = np.linspace(0, 6, 10000)
zb = 2 + 4*K**0.13 * np.exp(-2.5/K)
f = 1/np.sqrt(1 + (gamma/zb)**delta)
K_eff = K*f
zb_eff = 2 + 4*K_eff**0.13 * np.exp(-2.5/K_eff)

Delta = np.log(zb_eff/zb)-0.5*np.log(1+(gamma/zb_eff)**delta)



zeq = ((p*gamma**(delta/2))/K)**(2/(2+delta))
zb = 2+ 4*K**0.13 * np.exp(-2.5/K)
plt.plot(delta, Delta, 'r-', label='$\Delta$', linewidth=0.5)
#  plt.legend()
plt.xlabel('$\delta$')
plt.ylabel('$\Delta$')
plt.xlim([0, 6])
plt.ylim([-3, 1])
plt.show()