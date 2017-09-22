"""
This module defines a little utility class that can be used to derive
properties of HII regions from other properties.
"""

# Special check for readthedocs
import os
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

import numpy as np
import warnings
if not on_rtd:
    from scipy.optimize import brentq
else:
    def brentq(dummy1, dummy2, dummy3):
        pass

######################################################################
# Set some constants; change to cgs units
######################################################################

# Units and constants
try:
    from scipy.constants import c, m_e, m_p
    from scipy.constants import k as kB
    from scipy.constants import physical_constants as physcons
    c = c*1e2
    kB = kB*1e7
    m_e = m_e * 1e3
    m_p = m_p * 1e3
    eps0 = physcons['Rydberg constant times hc in J'][0] * 1e7
except:
    # This exception is to deal with readthedocs not having scipy
    c = 3.0e10
    kB = 1.38e-16
    m_e = 9.11e-28
    m_p = 1.67e-24
    
mH = m_e + m_p       # Hydrogen atom mass
alphaB = 2.59e-13    # Case B recombination coefficient
muH = 1.4            # Mean mass per H nucleus for standard cosmic composition
fe = 1.1             # Electrons per H nuclues
phi = 0.73           # Krumholz & Matzner dust phi parameter
TII = 1.0e4          # Fiducial temperature
psi = 3.2            # Krumholz & Matzner psi parameter
ftrap = 2.0          # Krumholz & Matzner trapping factor
Myr = 365.25*24.*3600.*1e6    # 1 Myr in seconds


######################################################################
# Define a little utility class
######################################################################
class windp(object):
    def __init__(self, fac=None):
        if fac is not None:
            self.fac = fac
    
    def windfac(self, omega):
        if omega < 1.0e6:
            return (1.0+omega)**(4./3.) - omega**(1./3.)*(4./3.+omega)
        else:
            # Series representation to avoid producing negative values
            return 2./9.*omega**(-2./3.)-4./81.*omega**(-5./3.) + \
                5./243.*omega**(-8./3.)

    def windfac1(self, omega):
        return self.windfac(omega)*(1.0+omega)**(1./6.)

    def windfac2(self, omega):
        return self.windfac(omega)*(9./4.)*omega**(2./3.)
        
    def windresid(self, logomega):
        return np.log(self.windfac(10.**logomega)/self.fac)

    def windresid1(self, logomega):
        return np.log(self.windfac1(10.**logomega)/self.fac)

    def windresid2(self, logomega):
        return np.log(self.windfac2(10.**logomega)/self.fac)

######################################################################
# This class takes inputs and derives nebular properties
######################################################################
class hiiregparam(object):
    
    ##################################################################
    def __init__(self, qH0, nII=None, r0=None, r1=None, U=None,
                 U0=None, Omega=None, t=None, n0=None, warn=True,
                 fix_quantity=None, r0safety=0.01):
        """
        Create an object to derive HII region parameters from inputs

        Parameters
           qH0 : float
              ionizing luminosity of stars, in photons / sec
           nII : float
              number density of H nuclei, in cm^-3
           rS : float
              Stromgren radius in cm^-3
           r0 : float
              inner radius in cm
           r1 : float
              outer radius in cm
           U : float
              volume-averaged ionization parameter
           U0 : float
              ionization parameter at inner radius
           Omega : float
              wind parameter
           t : float
              HII region age in Myr; if set, n0 and Omega must be set
              also
           n0 : float
              number density of H nuclei into which HII region is
              expanding, in cm^-3; if set, t and Omega must be set
              also
           warn : bool
              if True, warn and correct if parameter values are set
              outside physically possible range; if False, produce a
              ValueError when this happens
           fix_quantity : "nII" | "r0" | "r1" | "U" | "U0" | "Omega"
              if warn is True, this sepcifies which quantity is to be
              fixed in order to bring parameters into the allowed range
           r0safety : float
              if non-zero, when correcting quantities a fractional
              margin of safety of the specified value will be added to
              avoid having r0 = 0 exactly

        Notes
           Upon instantiation, the user must set either (1) exactly
           two of the six parameters (nII, r0, r1, U, U0, Omega); all
           pairwise combinations are allowed except (r0, U) and (r1,
           U0) -- these combinations do not define unique solutions;
           or (2) set the combination (t, n0, Omega). See slugpy
           documentation for an explanation of the underlying physics.
        """
        
        # Store properties
        self.qH0 = qH0
        self.nII_ = nII
        self.r0_ = r0
        self.r1_ = r1
        self.U_ = U
        self.U0_ = U0
        self.Omega_ = Omega
        self.warn = warn
        self.fix_quantity = fix_quantity
        self.r0safety = r0safety

        # If t and n0 are set, determine r1
        if t is not None and n0 is not None:
            self.r1_ = self.rKM(t, n0)

        # Make sure we have valid inputs
        self.check_params()

    ##################################################################
    def check_params(self):

        # Check we have the right number of parameters
        nparam = (self.nII_ is not None) + (self.r0_ is not None) + \
                 (self.r1_ is not None) + (self.U_ is not None) + \
                 (self.Omega_ is not None) + (self.U0_ is not None)
        if nparam != 2:
            raise ValueError(
                "hiiregparam: need exactly 2 of the " +
                "following: nII, r0, r1, U, U0 Omega")

        # Check for disallowed combinations (disallowed because the
        # solution is non-unique)
        if (self.U_ is not None) and (self.r0_ is not None):
            raise ValueError(
                "hiiregparam: combination U and r0 " +
                "not allowed because solution is non-unique")
        if (self.U0_ is not None) and (self.r1_ is not None):
            raise ValueError(
                "hiiregparam: combination U0 and r1 " +
                "not allowed because solution is non-unique")

        # Check for invalid numerical values of combinations and fix
        # if possible; if we are at a limit such that Omega = 1
        # exactly, flag it so that we don't generate errors later by
        # trying to solve for Omega numerically
        self.Omega_lim_ = False

        # Inner and outer radius
        if self.r0_ is not None and self.r1_ is not None:
            if self.r0_ >= self.r1_:
                if self.warn:
                    if self.fix_quantity != 'r0':
                        r1new = 1.0001*self.r0_
                        warnings.warn(
                            "hiiregparam: r1 must be >= r0; increasing " +
                            "r1 from {:e}".format(self.r1_) +
                            " to {:e}".format(r1new))
                        self.r1_ = r1new
                    else:
                        r0new = 0.9999*self.r0_
                        warnings.warn(
                            "hiiregparam: r1 must be >= r0; decreasing " +
                            "r0 from {:e}".format(self.r0_) +
                            " to {:e}".format(r0new))
                        self.r0_ = r0new
                else:
                    raise ValueError(
                        "hiiregparam: r1 must be >= r0")

        # Density and ionization parameter
        if self.nII_ is not None and self.U_ is not None:
            fac = self.U_**3*256*np.pi*c**3*fe / \
                  (81*alphaB**2*self.nII_*self.qH0)
            if fac == 1.0 and self.r0safety == 0.0:
                self.Omega_lim_ = True
            elif fac > 1 or (fac == 1.0 and self.r0safety > 0.0):
                Ulim = (81*alphaB**2*self.nII_*self.qH0 /
                           (256*np.pi*c**3*fe))**(1./3.)
                if self.warn:
                    if self.fix_quantity != "nII":
                        if self.r0safety > 0.0:
                            Unew = (1.0 - self.r0safety)*Ulim
                        else:
                            Unew = Ulim
                            self.Omega_lim_ = True
                        warnings.warn(
                            "hiiregparam: U too large for input "+
                            "value of nII; lowering from "+
                            "{:e}".format(self.U_)+" to "+
                            "{:e}".format(Unew))
                        self.U_ = Unew
                    else:
                        nIIlim =  self.U_**3*256*np.pi*c**3*fe / \
                                  (81*alphaB**2*self.qH0)
                        if self.r0safety > 0.0:
                            nIInew = (1.0+self.r0safety)*nIIlim
                        else:
                            nIInew = nIIlim
                            self.Omega_lim_ = True
                        warnings.warn(
                            "hiiregparam: nII too small for input "+
                            "value of U; increasing from {:e}".
                            format(self.nII_)+" to {:e}".
                            format(nIInew))
                        self.nII_ = nIInew
                else:
                    raise ValueError(
                        "hiiregparam: U too large for input "+
                        "value of nII; for nII = "+
                        "{:e}".format(self.nII_) +
                        ", maximum U is "+
                        "{:e}".format(Ulim))

        # Density and outer radius
        if self.r1_ is not None and self.nII_ is not None:
            rs = (3.*self.qH0/
                  (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
            if self.r1_ == rs and self.r0safety == 0.0:
                self.Omega_lim_ = True
            elif self.r1_ < rs or (self.r1_ == rs and self.r0safety > 0.0):
                nIIlim \
                    = np.sqrt(3.*self.qH0/
                              (4.*np.pi*self.r1_**3*alphaB*fe))
                if self.warn:
                    if self.fix_quantity != 'r1':
                        if self.r0safety > 0.0:
                            nIInew = (1.0+self.r0safety)*nIIlim
                        else:
                            nIInew = nIIlim
                            self.Omega_lim_ = True
                        warnings.warn(
                            "hiiregparam: nII too small for input "+
                            "value of r1; raising from "+
                            "{:e}".format(self.nII_)+" to "+
                            "{:e}".format(nIInew))
                        self.nII_ = nIInew
                    else:
                        r1lim = rs
                        if self.r0safety > 0.0:
                            r1new = (1.0+self.r0safety)*nIIlim
                        else:
                            r1new = r1lim
                            self.Omega_lim_ = True
                        warnings.warn(
                            "hiiregparam: r1 too small for input "+
                            "value of nII; raising from "+
                            "{:e}".format(self.r1_)+" to "+
                            "{:e}".format(r1new))
                        self.r1_ = r1new
                else:
                    raise ValueError(
                       "hiiregparam: nII too small for input "+
                        "value of r1; for r1 = " +
                        "{:e}".format(self.r1_) +
                        ", minimum nII is "+ 
                        "{:e}".format(nIIlim))

        # Outer radius and ionization parameter
        if self.r1_ is not None and self.U_ is not None:
            fac = 64.0*np.pi*c**2*fe*self.r1_*self.U_**2 / \
                  (81.*alphaB*self.qH0)
            if fac == 1 and self.r0safety == 0.0:
                self.Omega_lim_ = True
            elif fac > 1 or (fac == 1 and self.r0safety > 0.0):
                Ulim = np.sqrt(81.*alphaB*self.qH0 /
                               (64.*np.pi*c**2*fe*self.r1_))
                if self.warn:
                    if self.fix_quantity != 'r1':
                        if self.r0safety > 0.0:
                            Unew = (1.0-self.r0safety)*Ulim
                        else:
                            Unew = Ulim
                            self.Omega_lim_ = True
                        warnings.warn(
                            "hiiregparam: U too large for input "+
                            "value of r1; lowering from "+
                            "{:e}".format(self.U_) +
                            " to " +
                            "{:e}".format(Unew))
                        self.U_ = Unew
                    else:
                        r1lim = 81.*alphaB*self.qH0 / \
                                (64.0*np.pi*c**2*fe*self.U_**2)
                        if self.r0safety > 0.0:
                            r1new = (1.0 - self.r0safety)*r1lim
                        else:
                            r1new = r1lim
                            self.Omega_lim_ = True
                        warnings.warn(
                            "hiiregparam: r1 too large for input "+
                            "value of U; lowering from "+
                            "{:e}".format(self.r1_) +
                            " to " +
                            "{:e}".format(r1new))
                        self.r1_ = r1new
                else:
                    raise ValueError(
                        "hiiregparam: U too large for input "+
                        "value of r1; for r1 = " +
                        "{:e}".format(self.r1_) +
                        ", maximum U is " +
                        "{:e}".format(Ulim))

        # Inner and average ionization parameter
        if self.U_ is not None and self.U0_ is not None:
            if self.U0_ < 2.0*self.U_:
                Ulim = self.U0_*0.49822  # Corresponds to Omega = 10^6
                if self.warn:
                    warnings.warn(
                        "hiiregparam: U too large for input "+
                        " U0; maximum value of U is U0/2; lowering "+
                        " U to {:e}".format(Ulim))
                    self.U_ = Ulim
                else:
                    raise ValueError(
                        "hiiregparam: maximum value of U is U0/2")

    ##################################################################
    # Properties to derive all quantities
    ##################################################################

    ##################################################################
    @property
    def nII(self):
        if self.nII_ is not None:
            return self.nII_
        else:
            if self.r0_ is not None and self.r1_ is not None:
                 return np.sqrt(
                    3.*self.qH0/
                    (4.*np.pi*(self.r1_**3-self.r0_**3)*alphaB*fe))
            elif self.r0_ is not None and self.Omega_ is not None:
                rs = self.r0_ / self.Omega_**(1./3.)
                return np.sqrt(3.*self.qH0/
                               (4.*np.pi*alphaB*fe*rs**3))
            elif self.r0_ is not None and self.U0_ is not None:
                return self.qH0 / (4.*np.pi*self.r0**2*fe*c)
            elif self.r1_ is not None and self.U_ is not None:
                fac = 64.0*np.pi*c**2*fe*self.r1_*self.U_**2 / \
                      (81.*alphaB*self.qH0)
                wp = windp(fac=fac)
                if not self.Omega_lim_:
                    try:
                        Omega = 10.**brentq(wp.windresid1, -20, 20)
                    except ValueError:
                        raise ValueError(
                            "hiiregparam: failed to find solution "+
                            "for Omega at (r1 = {:e}, U = {:e}, qH0 = {:e})".
                            format(self.r1_, self.U_, self.qH0))
                else:
                    Omega = 0.0
                return 256.*np.pi*c**3*fe / \
                    (81.*alphaB**2*self.qH0) * \
                    (self.U_/wp.windfac(Omega))**3
            elif self.r1_ is not None and self.Omega_ is not None:
                rs = self.r1_ / (1+self.Omega_)**(1./3.)
                return np.sqrt(3.*self.qH0/
                               (4.*np.pi*alphaB*fe*rs**3))
            elif self.U_ is not None and self.U0_ is not None:
                wp = windp(fac=self.U/self.U0)
                try:
                    Omega = 10.**brentq(wp.windresid2, -20, 20)
                except ValueError:
                    raise ValueError(
                        "hiiregparam: failed to find solution "+
                        "for Omega at (U0 = {:e}, U = {:e})".
                        format(self.U0_, self.U_))
                fac = wp.windfac(Omega)
                return 256*fe*np.pi*c**3*(self.U_/fac)**3 / \
                    (81.*alphaB**2*self.qH0)
            elif self.U_ is not None and self.Omega_ is not None:
                wp = windp()
                fac = wp.windfac(self.Omega_)
                return 256*fe*np.pi*c**3*(self.U_/fac)**3 / \
                    (81.*alphaB**2*self.qH0)
            elif self.U0_ is not None and self.Omega_ is not None:
                return 36.*np.pi*c**3*fe*self.Omega_**2*self.U0_**3/ \
                    (alphaB*self.qH0)

    @nII.setter
    def nII(self, val):
        if self.nII_ is None:
            raise ValueError(
                "hiiregparam: cannot add a parameter after instantiation")
        else:
            self.nII_ = val
            self.check_params()


    ##################################################################
    @property
    def r0(self):
        if self.r0_ is not None:
            return self.r0_
        else:
            if self.nII_ is not None and self.r1_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
                return (self.r1_**3 - rs**3)**(1./3.)
            elif self.nII_ is not None and self.U_ is not None:
                fac = (self.U_**3*256*np.pi*c**3*fe / \
                       (81*alphaB**2*self.nII_*self.qH0))**(1./3.)
                wp = windp(fac=fac)
                if not self.Omega_lim_:
                    try:
                        Omega = 10.**brentq(wp.windresid, -20, 20)
                        rs = (3.*self.qH0/
                              (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
                        return Omega**(1./3.)*rs
                    except ValueError:
                        raise ValueError(
                            "hiiregparam: failed to find solution "+
                            "for Omega at (nII = {:e}, U = {:e}, qH0 = {:e})".
                            format(self.nII_, self.U_, self.qH0))
                else:
                    return 0.0
            elif self.nII_ is not None and self.U0_ is not None:
                return np.sqrt(self.qH0/
                               (4.*np.pi*fe*self.nII_*c*self.U0_))
            elif self.nII_ is not None and self.Omega_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
                return Omega**(1./3.)*rs
            elif self.r1_ is not None and self.U_ is not None:
                fac = 64.0*np.pi*c**2*fe*self.r1_*self.U_**2 / \
                      (81.*alphaB*self.qH0)
                wp = windp(fac=fac)
                if not self.Omega_lim_:
                    try:
                        Omega = 10.**brentq(wp.windresid1, -20, 20)
                    except ValueError:
                        raise ValueError(
                            "hiiregparam: failed to find solution "+
                            "for Omega at (r1 = {:e}, U = {:e}, qH0 = {:e})".
                            format(self.r1_, self.U_, self.qH0))
                else:
                    Omega = 0.0
                return self.r1_ * (Omega/(1.0+Omega))**(1./3.)
            elif self.r1_ is not None and self.Omega_ is not None:
                return (self.Omega_/(1.0+self.Omega_))**(1./3.) \
                    * self.r1_
            elif self.U_ is not None and self.U0_ is not None:
                nII = self.nII
                return np.sqrt(self.qH0/
                               (4.*np.pi*fe*nII*c*self.U0_))
            elif self.U_ is not None and self.Omega_ is not None:
                wp = windp()
                fac = wp.windfac(self.Omega_)
                rs = (fac/self.U_)**2*81*alphaB*self.qH0 / \
                     (64*np.pi*c**2*fe)
                return self.Omega_**(1./3.)*rs
            elif self.U0_ is not None and self.Omega_ is not None:
                wp = windp()
                fac = wp.windfac(self.Omega_)
                U = self.U0_*(9./4.)*self.Omega_**(2./3.)*fac
                rs = (fac/U)**2*81*alphaB*self.qH0 / \
                     (64*np.pi*c**2*fe)
                return self.Omega_**(1./3.)*rs

            
    @r0.setter
    def r0(self, val):
        if self.r0_ is None:
            raise ValueError(
                "hiiregparam: cannot add a parameter after instantiation")
        else:
            self.r0_ = val
            self.check_params()
            

    ##################################################################
    @property
    def r1(self):
        if self.r1_ is not None:
            return self.r1_
        else:
            if self.nII_ is not None and self.r0_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII**2))**(1./3.)
                return (rs**3+self.r0**3)**(1./3.)
            elif self.nII_ is not None and self.U_ is not None:
                fac = (self.U_**3*256*np.pi*c**3*fe / \
                       (81*alphaB**2*self.nII_*self.qH0))**(1./3.)
                wp = windp(fac=fac)
                if not self.Omega_lim_:
                    try:
                        Omega = 10.**brentq(wp.windresid, -20, 20)
                    except ValueError:
                        raise ValueError(
                            "hiiregparam: failed to find solution "+
                            "for Omega at (nII = {:e}, U = {:e}, qH0 = {:e})".
                            format(self.nII_, self.U_, self.qH0))
                else:
                    Omega = 0.0
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII**2))**(1./3.)
                return (1+Omega)**(1./3.)*rs
            elif self.nII_ is not None and self.U0_ is not None:
                Omega = np.sqrt(alphaB**2*self.nII_*self.qH0/
                                (36.*np.pi*c**3*fe)) * self.U0_**-1.5
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII**2))**(1./3.)
                return (1+Omega)**(1./3.)*rs                
            elif self.nII_ is not None and self.Omega_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII**2))**(1./3.)
                return (1+self.Omega_)**(1./3.)*rs
            elif self.r0_ is not None and self.U0_ is not None:
                nII = self.qH0/(4.*np.pi*self.r0_**2*fe*c*self.U0_)
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII**2))**(1./3.)
                return (rs**3+self.r0**3)**(1./3.)
            elif self.r0_ is not None and self.Omega_ is not None:
                rs = self.r0_ / self.Omega_**(1./3.)
                return (rs**3+self.r0**3)**(1./3.)
            elif self.U_ is not None and self.U0_ is not None:
                wp = windp(fac=self.U/self.U0)
                try:
                    Omega = 10.**brentq(wp.windresid2, -20, 20)
                except ValueError:
                    raise ValueError(
                        "hiiregparam: failed to find solution "+
                        "for Omega at (U0 = {:e}, U = {:e})".
                        format(self.U0_, self.U_))
                nII = 36.*np.pi*c**3*fe*Omega**2*self.U0_**3/ \
                      (alphaB*self.qH0)
                fac = wp.windfac(Omega)
                rs = (fac/self.U_)**2*81*alphaB*self.qH0 / \
                     (64*np.pi*c**2*fe)
                return (1.+Omega)**(1./3.)*rs
            elif self.U_ is not None and self.Omega_ is not None:
                wp = windp()
                fac = wp.windfac(self.Omega_)
                rs = (fac/self.U_)**2*81*alphaB*self.qH0 / \
                     (64*np.pi*c**2*fe)
                return (1.0+self.Omega_)**(1./3.)*rs
            elif self.U0_ is not None and self.Omega_ is not None:
                wp = windp()
                fac = wp.windfac(self.Omega_)
                U = self.U0*(9./4.)*self.Omega_**(2./3.)*fac
                rs = (fac/U)**2*81*alphaB*self.qH0 / \
                     (64*np.pi*c**2*fe)
                return (1.+self.Omega_)**(1./3.)*rs

    @r1.setter
    def r1(self, val):
        if self.r1_ is None:
            raise ValueError(
                "hiiregparam: cannot add a parameter after instantiation")
        else:
            self.r1_ = val

    ##################################################################
    @property
    def U(self):
        if self.U_ is not None:
            return self.U_
        else:
            if self.nII_ is not None and self.r0_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
                Omega = (self.r0_/rs)**3
                wp = windp()
                return (81.*alphaB**2*self.nII_*self.qH0 /
                        (256.*np.pi*c**3*fe))**(1./3.) \
                        * wp.windfac(Omega)
            elif self.nII_ is not None and self.r1_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
                Omega = (self.r1_/rs)**3 - 1.0
                wp = windp()
                return (81.*alphaB*self.qH0 /
                        (64.*np.pi*c**2*fe*rs))**(1./2.) \
                        * wp.windfac(Omega)
            elif self.nII_ is not None and self.U0_ is not None:
                Omega = np.sqrt(alphaB**2*self.nII_*self.qH0/
                                (36.*np.pi*c**3*fe)) * self.U0_**-1.5
                wp = windp()
                return (81.*alphaB**2*self.nII_*self.qH0 /
                        (256.*np.pi*c**3*fe))**(1./3.) \
                        * wp.windfac(Omega)                
            elif self.nII_ is not None and self.Omega_ is not None:
                wp = windp()
                return (81.*alphaB**2*self.nII_*self.qH0 /
                        (256.*np.pi*c**3*fe))**(1./3.) \
                        * wp.windfac(self.Omega_)
            elif self.r0_ is not None and self.r1_ is not None:
                nII = self.nII
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*nII**2))**(1./3.)
                Omega = (self.r0_/rs)**3
                wp = windp()
                return (81.*alphaB**2*nII*self.qH0 /
                        (256.*np.pi*c**3*fe))**(1./3.) \
                        * wp.windfac(Omega)
            elif self.r0_ is not None and self.U0_ is not None:
                nII = self.qH0/(4.*np.pi*self.r0_**2*fe*c*self.U0_)
                Omega = np.sqrt(alphaB**2*nII*self.qH0/
                                (36.*np.pi*c**3*fe)) * self.U0_**-1.5
                wp = windp()
                return (81.*alphaB**2*nII*self.qH0 /
                        (256.*np.pi*c**3*fe))**(1./3.) \
                        * wp.windfac(Omega)
            elif self.r0_ is not None and self.Omega_ is not None:
                rs = self.r0_/self.Omega_**(1./3.)
                wp = windp()
                return (81.*alphaB*self.qH0 /
                        (64.*np.pi*c**2*fe*rs))**(1./2.) \
                        * wp.windfac(self.Omega_)
            elif self.r1_ is not None and self.Omega_ is not None:
                rs = self.r1_/(1.0+self.Omega_)**(1./3.)
                wp = windp()
                return (81.*alphaB*self.qH0 /
                        (64.*np.pi*c**2*fe*rs))**(1./3.) \
                        * wp.windfac(self.Omega_)
            elif self.U0_ is not None and self.Omega_ is not None:
                wp = windp()
                fac = wp.windfac(self.Omega_)
                return self.U0*(9./4.)*self.Omega_**(2./3.)*fac

    @U.setter
    def U(self, val):
        if self.U_ is None:
            raise ValueError(
                "hiiregparam: cannot add a parameter after instantiation")
        else:
            self.U_ = val
            self.check_params()

    ##################################################################
    @property
    def Omega(self):
        if self.Omega_ is not None:
            return self.Omega_
        else:
            if self.nII_ is not None and self.r0_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
                return (self.r0_/rs)**(1./3.)
            elif self.nII_ is not None and self.r1_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
                return (self.r1_/rs)**(1./3.) - 1.0
            elif self.nII_ is not None and self.U_ is not None:
                fac = (self.U_**3*256*np.pi*c**3*fe / \
                      (81*alphaB**2*self.nII_*self.qH0))**(1./3.)
                wp = windp(fac=fac)
                if not self.Omega_lim_:
                    try:
                        Omega = 10.**brentq(wp.windresid, -20, 20)
                    except ValueError:
                        raise ValueError(
                            "hiiregparam: failed to find solution "+
                            "for Omega at (nII = {:e}, U = {:e}, qH0 = {:e})".
                            format(self.nII_, self.U_, self.qH0))
                else:
                    Omega = 0.0
                return Omega
            elif self.nII_ is not None and self.U0_ is not None:
                return np.sqrt(alphaB**2*self.nII_*self.qH0/
                               (36.*np.pi*c**3*fe)) * self.U0_**-1.5
            elif self.r0_ is not None and self.r1_ is not None:
                return self.r0_**3 / (self.r1_**3 - self.r0_**3)
            elif self.r0_ is not None and self.U0_ is not None:
                nII = self.qH0/(4.*np.pi*self.r0_**2*fe*c*self.U0_)
                return np.sqrt(alphaB**2*nII*self.qH0/
                               (36.*np.pi*c**3*fe)) * self.U0_**-1.5
            elif self.r1_ is not None and self.U_ is not None:
                fac = 64.0*np.pi*c**2*fe*self.r1_*self.U_**2 / \
                      (81.*alphaB*self.qH0)
                wp = windp(fac=fac)
                if not self.Omega_lim_:
                    try:
                        Omega = 10.**brentq(wp.windresid1, -20, 20)
                    except ValueError:
                        raise ValueError(
                            "hiiregparam: failed to find solution "+
                            "for Omega at (r1 = {:e}, U = {:e}, qH0 = {:e})".
                            format(self.r1_, self.U_, self.qH0))
                else:
                    Omega = 0.0
                return Omega
            elif self.U_ is not None and self.U0_ is not None:
                wp = windp(fac=self.U_/self.U0_)
                try:
                    Omega = 10.**brentq(wp.windresid2, -20, 20)
                    return Omega
                except ValueError:
                    raise ValueError(
                        "hiiregparam: failed to find solution "+
                        "for Omega at (U0 = {:e}, U = {:e})".
                        format(self.U0_, self.U_))
            
    @Omega.setter
    def Omega(self, val):
        if self.Omega_ is None:
            raise ValueError(
                "hiiregparam: cannot add a parameter after instantiation")
        else:
            self.Omega_ = val
            self.check_params()

    ##################################################################
    @property
    def U0(self):
        if self.U0_ is not None:
            return self.U0_
        elif self.Omega_lim_:
            return np.inf
        else:
            if self.nII_ is not None and self.r0_ is not None:
                return self.qH0 / \
                    (4.*np.pi*self.r0_**2*fe*self.nII_*c)
            elif self.nII_ is not None and self.r1_ is not None:
                rs = (3.*self.qH0/
                      (4.*np.pi*alphaB*fe*self.nII_**2))**(1./3.)
                Omega = (self.r1_/rs)**3 - 1.0
                return (alphaB**2*self.nII_*self.qH0 /
                        (36.*np.pi*c**3*fe))**(1./3.) / Omega**(2./3.)
            elif self.nII_ is not None and self.U_ is not None:
                fac = (self.U_**3*256*np.pi*c**3*fe / \
                      (81*alphaB**2*self.nII_*self.qH0))**(1./3.)
                wp = windp(fac=fac)
                try:
                    Omega = 10.**brentq(wp.windresid, -20, 20)
                except ValueError:
                    raise ValueError(
                        "hiiregparam: failed to find solution "+
                        "for Omega at (nII = {:e}, U = {:e}, qH0 = {:e})".
                        format(self.nII_, self.U_, self.qH0))
                return (alphaB**2*self.nII_*self.qH0 /
                        (36.*np.pi*c**3*fe))**(1./3.) / Omega**(2./3.)
            elif self.nII_ is not None and self.Omega_ is not None:
               return (alphaB**2*self.nII_*self.qH0 /
                        (36.*np.pi*c**3*fe))**(1./3.) / \
                        self.Omega_**(2./3.)
            elif self.r0_ is not None and self.r1_ is not None:
                nII = self.nII
                return self.qH0 / \
                    (4.*np.pi*self.r0_**2*fe*nII*c)
            elif self.r0_ is not None and self.Omega_ is not None:
                U = self.U
                wp = windp()
                fac = wp.windfac2(self.Omega_)
                return U/fac
            elif self.r1_ is not None and self.U_ is not None:
                Omega = self.Omega
                wp = windp()
                fac = wp.windfac2(Omega)
                return self.U_ / fac
            elif self.r1_ is not None and self.Omega_ is not None:
                U = self.U
                wp = windp()
                fac = wp.windfac2(self.Omega_)
                return U/fac
            elif self.U_ is not None and self.Omega_ is not None:
                wp = windp()
                fac = wp.windfac2(self.Omega_)
                return self.U_/fac
            
    @U0.setter
    def U0(self, val):
        if self.U0_ is None:
            raise ValueError(
                "hiiregparam: cannot add a parameter after instantiation")
        else:
            self.U0_ = val
            self.check_params()

    ##################################################################
    def rS(self):
        """
        Return the Stromgren radius 

        Parameters
           None

        Returns
           rS : float
              The value of rS in cm
        """
        
        rS = (3.0*self.qH0/(4.0*np.pi*alphaB*fe*self.nII**2))**(1./3.)
        return rS
   
    ##################################################################
    def rch(self):
        """
        Return the Krumholz & Matzner characteristic radius, which
        defines the radius inside which radiation pressure is
        dynamically significant. 

        Parameters
           None

        Returns
           rch : float
              The value of r_ch in cm
        """
        
        rch = alphaB/(12.*np.pi*phi) * (eps0/(2.0*fe*kB*TII))**2 \
              * ftrap**2 * psi**2 * self.qH0/ c**2
        return rch

    ##################################################################
    def zeta(self):
        """
        Return the Krumholz & Matzner zeta parameter, defined as
        r_ch/r_1. Values of zeta >> 1 imply that radiation pressure is
        dynamically significant.

        Parameters
           None

        Returns
           zeta : float
              The value of zeta for the current HII region
        """
        
        rch = alphaB/(12.*np.pi*phi) * (eps0/(2.0*fe*kB*TII))**2 \
              * ftrap**2 * psi**2 * self.qH0/ c**2
        return rch / self.r1

    ##################################################################
    def rKM(self, t, n):
        """
        Compute the radius predicted by Krumholz & Matzner (2009) for
        a specified age and ambient density

        Parameters
           t : float
              HII region age in Myr
           n : float
              density of region into which HII region is expanding, in
              cm^-3

        Returns
           r : float
              radius predicted by Krumholz & Matnzer (2009)
        """
        rch = self.rch()
        tch = np.sqrt(
            4.*np.pi*muH*mH*n*c*rch**4 /
            (3.*ftrap*self.qH0*psi*eps0))
        tau = t*Myr/tch
        xrad = (2.*tau**2)**0.25
        xgas = (49.*tau**2/36.)**(2./7.)
        return rch * (xrad**3.5 + xgas**3.5)**(2./7.)
