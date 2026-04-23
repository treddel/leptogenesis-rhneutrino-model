import numpy as np
import matplotlib.pyplot as plt

K = 100
gamma = 50
p = 10

delta = np.linspace(0, 6, 10000)

zeq = ((p*gamma**(delta/2))/K)**(2/(2+delta))
zb = 2+ 4*K**0.13 * np.exp(-2.5/K)
plt.plot(delta, zeq, 'r-', label='D', linewidth=0.5)
plt.xlabel('$\delta$')
plt.ylabel('$z_{eq}$')
plt.xlim([0, 6])
plt.ylim([0, 10])
plt.show()