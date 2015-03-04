# w_scan.py
import sys, os

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here, '../src')))

from phasematching import *
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
from matplotlib import pyplot as plt


num = 100
pulse_scan = np.linspace(10e-15, 300e-15, num)

harm = np.array(np.zeros(num))
Lcoh = np.array(np.zeros(num))
eta = np.array(np.zeros(num))
Labs = np.array(np.zeros(num))
Lmed = np.array(np.zeros(num))

inten = .5
for i in range(num):
    sim = phase_matching('Xe', 17,  inten ,  pulse_scan[i], 17e-6, 1070e-9,  .17, 0.35e-3,  100, 0, 'on')
    harm[i], Lcoh[i], eta[i], Labs[i], Lmed[i] = sim.int_harmonic_yield()

fig, ax = plt.subplots(2,2, figsize = (15,8))
ax[0,0].semilogy(pulse_scan * 10 ** 15, harm, 'k-', linewidth = 2)
ax[0,0].set_xlabel('Pulse Duration [fs]')
ax[0,0].set_ylabel('Harmonic Yield [arb.]')
ax[0,0].set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax[0,0].grid()




ax[0,1].plot(pulse_scan * 10 ** 15, Lcoh * 10 ** 6, 'k-', linewidth = 2, label = 'Coherence Length')
ax[0,1].plot(pulse_scan * 10 ** 15, Lmed * 10 ** 6, 'r-', linewidth = 2, label = 'Medium Length')
ax[0,1].plot(pulse_scan * 10 ** 15, Labs * 10 ** 6, 'b-', linewidth = 2, label = 'Absorption Length')
ax[0,1].set_xlabel('Pulse Duration [fs]')
ax[0,1].set_ylabel('Length [um]')
ax[0,1].legend(prop={'size':12})
ax[0,1].grid()


ax[1,0].plot(pulse_scan * 10 ** 15, eta, 'm-', linewidth = 2, label = 'Ionization Fraction')
ax[1,0].set_xlabel('Pulse Duration [fs]')
ax[1,0].set_ylabel('Ionization Fraction')
ax[1,0].grid()


ax[1,1].plot(pulse_scan * 10 ** 15, Lcoh / Labs, 'k-', linewidth = 2, label = 'L_coh / L_abs ( > 3)')
ax[1,1].plot(pulse_scan * 10 ** 15, Lmed / Labs, 'r-', linewidth = 2, label = 'L_med / L_abs ( > 5)')
# ylim(.01,1)
ax[1,1].set_xlabel('Pulse Duration [fs]')
ax[1,1].set_ylabel('Ratio')

ax[1,1].legend(prop={'size':12})
ax[1,1].grid()
fig.set_tight_layout(True)


fig, ax = plt.subplots(figsize = (15,8))
ax.plot(pulse_scan * 10 ** 15, inten * 10 ** 15 * np.pi/2 * pulse_scan * 154e6 *(17e-6) ** 2, 'k-', linewidth = 2, label = ' Avg. Power Required')

ax.set_xlabel('Pulse Duration [fs]')
ax.set_ylabel('Power [kW]')
ax.set_title(sim.Atom + ' at ' + str(sim.q) + ' Harmonic')
ax.legend(prop={'size':12})
ax.grid()

fig.set_tight_layout(True)
plt.show()