import numpy as np
import matplotlib.pyplot as plt

N_arb = [769, 1153, 3889, 10369, 12289, 17497]
t_arb = [15799, 24940, 111777, 315829, 351998, 587744]

N_arb = np.array(N_arb)
t_arb = np.array(t_arb)

N_pw2 = [512, 1024, 8192, 16384]
t_pw2 = [3571, 7465, 84758, 185806]

N_pw2 = np.array(N_pw2)
t_pw2 = np.array(t_pw2)

# use log scale for both axes

plt.figure()
plt.title('FFT performance')
plt.xscale('log')
plt.yscale('log')

plt.plot(N_arb, t_arb, 'o-', label='Arbitrary')
plt.plot(N_pw2, t_pw2, 'o-', label='Power of 2')
plt.xlabel('Cyclotomic order')
plt.ylabel('Time (ns)')
plt.legend()
plt.show()