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


def model_ivp_2(z, N):
    N1 = N[0]
    Nb = N[1]
    K_100 = 10**2
    ep = 10**-6

    dN1dz_100 = -(K_100*z*(kn(1, z)/kn(2, z)))*(N1 - (0.375* z**2 *kn(2,z)))
    dNbdz_100 = -ep*(K_100*z*(kn(1, z)/kn(2, z)))*(N1 - (0.375* z**2 *kn(2,z))) - (0.25*K_100* z**3 *kn(3, z)*Nb)

    dNdz_100 = [dN1dz_100, dNbdz_100]
    return dNdz_100


def model_ivp_4(z, N):
    N1 = N[0]
    Nb = N[1]
    K_4 = 10**-4
    ep = 10**-6

    dN1dz_4 = -(K_4*z*(kn(1, z)/kn(2, z)))*(N1 - (0.375* z**2 *kn(2,z)))
    dNbdz_4 = -ep*(K_4*z*(kn(1, z)/kn(2, z)))*(N1 - (0.375* z**2 *kn(2,z))) - (0.25*K_4* z**3 *kn(3, z)*Nb)

    dNdz_4 = [dN1dz_4, dNbdz_4]
    return dNdz_4


N0 = [0.75, 0]
z_span = [10**-4, 10**3]
z = np.linspace(z_span[0], z_span[1], 1000)

N_ivp_2 = solve_ivp(model_ivp_2, z_span, N0, atol=1e-14, rtol=1e-10)
N_ivp_4 = solve_ivp(model_ivp_4, z_span, N0, atol=1e-14, rtol=1e-10)
N_ivp = solve_ivp(model_ivp, z_span, N0, atol=1e-14, rtol=1e-10)
N1eq = 0.375 * z**2 * kn(2, z)
kappa_2 = 4 / 3 * (N0[0] - N_ivp_2.y[0]) * 10 ** -2
kappa_4 = 4/3 * (N0[0] - N_ivp_4.y[0]) * 10**-2
kappa = 4/3 * (N0[0] - N_ivp.y[0]) * 10**-2

plt.plot(N_ivp.t, np.abs(kappa), 'b--', label='$|\kappa| x $10$^{-2}$', linewidth=0.5)
plt.plot(N_ivp.t, np.abs(N_ivp.y[0]), 'b-', label='$\eta_{N_1}$', linewidth=0.5)
plt.plot(z, N1eq, 'k--', label='$\eta_{N_1}^{eq}$', linewidth=0.5)
plt.plot(N_ivp_2.t, np.abs(kappa_2), 'r--', linewidth=0.5)
plt.plot(N_ivp_2.t, np.abs(N_ivp_2.y[0]), 'r-', linewidth=0.5)
plt.plot(N_ivp_4.t, np.abs(kappa_4), 'g--', linewidth=0.5)
plt.plot(N_ivp_4.t, np.abs(N_ivp_4.y[0]), 'g-', linewidth=0.5)

plt.legend()
plt.xlabel('z = M$_1$ / T')
plt.ylabel('')
plt.xscale('log')
plt.yscale('log')
plt.xlim([10**-1, 10**3])
plt.ylim([10**-4, 10**1])
plt.show()

plt.close('all')
