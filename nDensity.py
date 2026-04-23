import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt

def model_ivp(z, N):
    N1 = N[0]
    Nb = N[1]
    K = 10**-2
    ep = 10**-6

    dN1dz = -(K*z*(kn(1, z)/kn(2, z)))*(N1 - (0.375* z**2 *kn(2,z)))
    dNbdz = -ep*(K*z*(kn(1, z)/kn(2, z)))*(N1 - (0.375* z**2 *kn(2,z))) - (0.25*K* z**3 *kn(3, z)*Nb)

    dNdz = [dN1dz, dNbdz]
    return dNdz


N0 = [0, 0]
z_span = [10**-4, 10**2]
z = np.linspace(z_span[0], z_span[1], 1000)

N_ivp = solve_ivp(model_ivp, z_span, N0, atol=1e-14, rtol=1e-10)
N1eq = 0.375 * z**2 * kn(2, z)

# solve_ivp plotting
plt.plot(N_ivp.t, 'r-', N_ivp.y[0], label='$N_1$', linewidth=0.5)
plt.plot(N_ivp.t, 'r-', np.abs(N_ivp.y[1]), label='$N_{B-L}$', linewidth=0.5)
plt.plot(z, N1eq, 'k--', label='$N_1^{eq}$', linewidth=0.5)
# plt.legend()
plt.xlabel('z = T / m$_1$')
plt.ylabel('N')  # ($\eta_\gamma = 1$ normalisation)
# plt.title("Weak Washout Evolutions (K = 10$^{-2}$)")
plt.xscale('log')
plt.yscale('log')
plt.xlim([10**-2, 10**2])
plt.ylim([10**-11, 10**1])
plt.show()

plt.close('all')
