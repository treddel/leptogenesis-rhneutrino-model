import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt

K = 100
ep = 10 ** -6

y_span = [10**-11, 10**4]
z_span = [10**-4, 10**2]
z = np.linspace(z_span[0], z_span[1], 160764)
y_line = np.linspace(y_span[0], 0.841, 1000)
zb = 9.206
zb_line = np.zeros_like(y_line) + zb

Wid = 0.25* K * z**3*kn(1,z)
zb_cond = (kn(2,z)/kn(1,z)) - 3/z

plt.plot(z, Wid, label='$W_{ID}$')
plt.plot(z, zb_cond)
plt.plot(zb_line, y_line, '--')
plt.legend()
plt.xlabel('z = T / m$_1$')
plt.ylabel('')
#  plt.title("Washout ")
plt.xscale('log')
plt.yscale('log')
plt.xlim([10**-2, 10**2])
plt.ylim([10**-11, 10**4])
plt.show()
