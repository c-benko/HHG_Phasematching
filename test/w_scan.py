# w_scan.py
import sys, os

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here, '../src')))

from phasematching import *
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
from matplotlib import pyplot as plt

num = 100
wscan = np.linspace(10e-6, 30e-6, num)

harm = np.array(np.zeros(num))
Lcoh = np.array(np.zeros(num))
eta = np.array(np.zeros(num))
Labs = np.array(np.zeros(num))
Lmed = np.array(np.zeros(num))

inten = .5
pw = 120e-15
for i in range(num):
    sim = phase_matching('Xe', 17,  .6 ,  pw, wscan[i], 1070e-9,  .17, 0.15e-3,  100, 0, 'on')
    harm[i], Lcoh[i], eta[i], Labs[i], Lmed[i] = sim.int_harmonic_yield()

fig, ax = plt.subplots(2,2, figsize = (15,8))

ax[0,0].semilogy(wscan * 10 ** 6, harm, 'k-', linewidth = 2)
ax[0,0].set_xlabel('Spot Size [um]')
ax[0,0].set_ylabel('Harmonic Yield [arb.]')
ax[0,0].set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax[0,0].grid()




ax[0,1].plot(wscan * 10 ** 6, Lcoh * 10 ** 3, 'k-', linewidth = 2, label = 'Coherence Length')
ax[0,1].plot(wscan * 10 ** 6, Lmed * 10 ** 3, 'r-', linewidth = 2, label = 'Medium Length')
ax[0,1].plot(wscan * 10 ** 6, Labs * 10 ** 3, 'b-', linewidth = 2, label = 'Absorption Length')
ax[0,1].set_xlabel('Spot Size [um]')
ax[0,1].set_ylabel('Length [um]')
ax[0,1].set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax[0,1].legend(prop={'size':12})
ax[0,1].grid()

ax[1,0].plot(wscan * 10 ** 6, eta, 'm-', linewidth = 2, label = 'Ionization Fraction')
ax[1,0].set_xlabel('Spot Size [um]')
ax[1,0].set_ylabel('Ionization Fraction')
ax[1,0].set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax[1,0].grid()


ax[1,1].plot(wscan * 10 ** 6, Lcoh / Labs, 'k-', linewidth = 2, label = 'L_coh / L_abs ( > 3)')
ax[1,1].plot(wscan * 10 ** 6, Lmed / Labs, 'r-', linewidth = 2, label = 'L_med / L_abs ( > 5)')
ax[1,1].set_xlabel('Spot Size [um]')
ax[1,1].set_ylabel('Ratio')
ax[1,1].set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax[1,1].legend(prop={'size':12})
ax[1,1].grid()
fig.set_tight_layout(True)

fig, ax = plt.subplots()
ax.plot(wscan * 10 ** 6, inten * 10 ** 15 * np.pi/2 * pw * 154e6 * wscan ** 2, 'k-', linewidth = 2, label = ' Avg. Power Required')

ax.set_xlabel('Spot Size [um]')
ax.set_ylabel('Power [kW]')
ax.set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax.legend(prop={'size':12})
ax.grid()
fig.set_tight_layout(True)

plt.show()