# phasematching.py
# Absorption limited case.

# Need to add a part that calculates the nonlinear phase shift to make sure it never reaches pi.
# also add energy loss per Tom's simulation (Eq 5.9 of SM)

# one one to take into account the cavity is to use injected power as input, then from 
# tom's equations I can calculate the reduced intensity due to energy transfer to plasma
# and the nonlinear phase shift. This should clamp how much I can inject into a cavity.


# my imports
from atom import *
from gas import *
from laser import *
from ADK import *

class phase_matching:
    '''
    Contains a calculation of phase matching conditions in the 
    absorption-limited case.

    References:
    E. Constant et al. Phys. Rev. Lett.  82, 1668–1671 (1999).
    T. K. Allison. UC Berkeley Thesis (2010).
    S. Hadrich et al. Nat. Photon. 8, 779–783 (2014).

    '''
    def __init__(self, Atom = 'Xe', q = 17, Intensity = 0.5 , Pulse_FWHM = 120e-15, 
        Spot = 17e-6, Lam = 1070e-9, Pressure = .1, Length = 0.1e-3, Temperature = 100, Z = 0, SS = 'on'):

        ## input quantities
        self.Atom = Atom
        self.q = q
        self.Intensity = Intensity      
        self.Pulse_FWHM = Pulse_FWHM
        self.Spot = Spot 
        self.Lam = Lam
        self.Pressure = Pressure
        self.Length = Length 
        self.Temperature = Temperature
        self.Z = Z
        self.SS = SS

        ## generated quantities used in functions below
        # atom
        self.at = atom(self.Atom, self.Lam , self.Pressure , self.Temperature)
        # abs cross section
        self.cross = self.at.xuv_absorption(self.q * 1240 / (self.Lam * 10 ** 9))
        # Ip
        self.Ip = self.at.adk_params()['Ip']
        # laser
        self.las = laser(self.Pulse_FWHM, self.Intensity, self.Lam, self.Spot)
        # gas jet
        self.targ = gas(self.Pressure, self.Temperature, self.Length)
        # adk class
        self.adk = ADK(self.Atom, self.Intensity, self.Pulse_FWHM)

    def dipole(self, inten):
        '''
        Returns the dipole. 
        '''
         
        Up = 9.33 * inten * (self.Lam * 10 ** 6)**2 
        harm = self.q * 1240 / (self.Lam * 10 ** 9) 
        Icut =  (harm - self.Ip) / (3.14 * 9.33 * (self.Lam * 10 ** 6) ** 2)
        if inten < Icut:
            Aq = ( inten / Icut ) ** 10.6
        else:
            Aq = ( inten / Icut ) ** 4.6
        return Aq

    def abs_length(self):
        '''
        Returns the absorption length in [m].
        '''

        g = gas(self.Pressure, self.Temperature, self.Length)
        return 1 / ( self.cross *  g.density )

    def abs_cav(self):
        '''
        Returns the absorption length in [m] for an XUV beam propagating
        in a residual background pressure in the vacuum chamber.

        The propagation distance is assumed to be 1 M. 

        1 atm of backing pressure is approximately 1.3e-5 to 1.3e-6 atm
        of residual gas.
        '''

        g = gas(self.Pressure * 1.3e-6, 300, 1)
        return np.exp(-1 * self.cross *  g.density * 1 ) 


    def coh_length(self, eta):
        '''
        Returns the coherence length in [m].
        '''
        # constants
        re = 2.8179 * 10 ** -15 #classical electron radius
        kb = 1.3806488 * 10 ** -23 #Boltzmann constant

        # index difference
        dn, eta_crit = self.at.eta_crit(self.q * 1240 / (self.Lam * 10 ** 9))

        # rayleigh range
        zr = self.las.zr

        # wavevector mismatch
        dk_atomic = self.q * 2 * np.pi / self.Lam * dn
        dk_plasma = -1*self.q * self.targ.density * re * self.Lam * eta
        dk_gouy = -1*self.q / zr / (1 + self.Z **2 / zr ** 2 )
        dk_dip = 2 * self.at.alpha1 * self.Intensity / zr * (self.Z / zr) / (1 + self.Z **2 / zr ** 2 )
    
        return np.pi / abs(dk_atomic + dk_plasma + dk_gouy + dk_dip)

    def harmonic_yield(self):
        '''
        Returns the harmonic yield as a function of time. 
        This is required to perform the integrated signal. 

        '''

        # time span
        num = 100
        tspan = np.linspace(-2 * self.Pulse_FWHM, 2 * self.Pulse_FWHM, num)
        dt = tspan[2] - tspan[1]

        # parameters
        Labs = self.abs_length()
        Lmed = self.targ.length
        Dens = self.targ.density

        # bins for outputs
        wbar = np.array(np.zeros(num))
        harm_yield = np.array(np.zeros(num))
        eta = np.array(np.zeros(num))
        Lcoh = np.array(np.zeros(num))

        # Steady-state ionization fraction
        if self.SS == 'on':
            eta_ss = self.adk.steady_state(self.adk.ion_frac_adk_TL(), self.at.kp())
        else:
            eta_ss = self.adk.steady_state(self.adk.ion_frac_adk_TL(), 1e6)

        for i in range(num):
            wbar[i] = self.adk.wbar_adk_TL(self.las.pulse(tspan[i]))
            eta[i] =  1 - (1 - eta_ss) * np.exp(-1 * sum(wbar) * dt) 
            Lcoh[i] = self.coh_length(eta[i])

            fac1 = 4 * Labs ** 2 * Dens ** 2 * (1 - eta[i]) ** 2  * self.dipole( self.las.pulse(tspan[i]) )
            fac2 = 1 / (1 + 4 * np.pi ** 2 * Labs ** 2 / Lcoh[i] ** 2)
            fac3 = 1 + np.exp(-Lmed / Labs) - 2 * np.cos(np.pi * Lmed / Lcoh[i]) * np.exp(-Lmed / 2 / Labs)
            harm_yield[i] =  fac1 * fac2 * fac3

        return tspan, self.las.pulse(tspan), Lcoh, eta, harm_yield, num

    def int_harmonic_yield(self):
        '''
        Returns the integrated harmonic yield. 
        This is required to perform the integrated signal. 
        '''
        tspan, pulse, Lcoh, eta, harm_yield, num = self.harmonic_yield()
        dt = tspan[2] - tspan[1]
        return sum(harm_yield) * dt * self.abs_cav(), Lcoh[int(num/2)], eta[-1], self.abs_length(), self.targ.length
