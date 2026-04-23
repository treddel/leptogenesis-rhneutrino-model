import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt

K = 10**-2
ep = 10 ** -6

def model_ivp(z, N):
    N1 = N[0]
    Nb = N[1]
    N1eq = N[2]

    dN1dz = -(K*z*(kn(1, z)/kn(2, z)))*(N1 - (0.375* z**2 *kn(2,z)))
    dNbdz = -ep*(K*z*(kn(1, z)/kn(2, z)))*(N1 - (0.375* z**2 *kn(2,z))) - (0.25*K* z**3 *kn(3, z)*Nb)
    dN1eqdz = -0.375 * z ** 2 * kn(1, z)

    dNdz = [dN1dz, dNbdz, dN1eqdz]
    return dNdz


N0 = [0, 0, 0.75]
z_span = [10**-4, 10**2]
z = np.linspace(z_span[0], z_span[1], 160685)
zbar = np.zeros_like(z)
zb = 9.206  # derived in zb.py


for i in range(len(z)):
    zbar[i] = min(z[i], zb)


N_ivp = solve_ivp(model_ivp, z_span, N0, atol=1e-14, rtol=1e-10)
N1eq = 0.375 * z**2 * kn(2, z)
Wid = 0.25* K * z**3*kn(1,z)
Dz = K * z**2 * (kn(1,z)/kn(2,z))
delta = N_ivp.y[0] - N_ivp.y[2]
delta_approx = (1 + (K * z**3/(15/4)+2*z))**-1 * (3/16) * z**3 * kn(1,z)
kappa = (1 + ((K**2 * zbar * z**5)/75))**-1 * 2*K/75 * z**5
Nbl = 0.75*ep * kappa


plt.plot(z, Dz, 'k--', label='$Dz$', linewidth=0.5)  # Decay term
plt.plot(z, Wid, 'k--', label='$W_{ID}$', linewidth=0.5)  # Washout term
# plt.plot(z, delta_approx, label='$\Delta$')
plt.plot(N_ivp.t, delta, 'r-', label='$\Delta$', linewidth=0.5)  #N1 - N1eq
plt.plot(N_ivp.t, np.abs(N_ivp.y[1]), 'b-', label='$N_{B-L}$', linewidth=0.5)  # N_B-L
plt.legend()
plt.xlabel('z = M$_1$ / T')
plt.ylabel('')
# plt.title("Strong Washout Evolutions (K = 10$^2$)")
plt.xscale('log')
plt.yscale('log')
plt.xlim([10**-2, 10**2])
plt.ylim([10**-11, 10**4])
plt.show()

