# i_scan.py
import sys, os

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here, '../src')))

from phasematching import *
import matplotlib
from matplotlib import pyplot as plt

num = 100
pscan = np.linspace(0.001, 3, num)

harm = np.array(np.zeros(num))
Lcoh = np.array(np.zeros(num))
eta = np.array(np.zeros(num))
Labs = np.array(np.zeros(num))
Lmed = np.array(np.zeros(num))

for i in range(num):
    sim = phase_matching('Xe', 17,  0.5 ,  120e-15, 17e-6, 1070e-9,  pscan[i], 0.10e-3,  100, 0, 'on')
    harm[i], Lcoh[i], eta[i], Labs[i], Lmed[i] = sim.int_harmonic_yield()

fig, ax = plt.subplots(2,2, figsize = (15,8))
matplotlib.rcParams.update({'font.size': 20})

ax[0,0].plot(pscan, harm, 'k-', linewidth = 2)
# ylim(.01,1)
ax[0,0].set_xlabel('Pressure [atm.]')
ax[0,0].set_ylabel('Harmonic Yield [arb.]')
ax[0,0].set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax[0,0].grid()




ax[0,1].plot(pscan, Lcoh * 10 ** 3, 'k-', linewidth = 2, label = 'Coherence Length')
ax[0,1].plot(pscan, Lmed * 10 ** 3, 'r-', linewidth = 2, label = 'Medium Length')
ax[0,1].plot(pscan, Labs * 10 ** 3, 'b-', linewidth = 2, label = 'Absorption Length')
ax[0,1].set_xlabel('Pressure [atm.]')
ax[0,1].set_ylabel('Length [m]')
ax[0,1].legend(prop={'size':12})
ax[0,1].grid()



ax[1,0].plot(pscan, eta, 'm-', linewidth = 2, label = 'Ionization Fraction')

ax[1,0].set_xlabel('Pressure [atm.]')
ax[1,0].set_ylabel('Ionization Fraction')
ax[1,0].grid()



ax[1,1].plot(pscan, Lcoh / Labs, 'k-', linewidth = 2, label = 'L_coh / L_abs ( > 3)')
ax[1,1].plot(pscan, Lmed / Labs, 'r-', linewidth = 2, label = 'L_med / L_abs ( > 5)')
ax[1,1].set_xlabel('Pressure [atm.]')
ax[1,1].set_ylabel('Ratio')
ax[1,1].legend(prop={'size':12})
ax[1,1].grid()

plt.show()