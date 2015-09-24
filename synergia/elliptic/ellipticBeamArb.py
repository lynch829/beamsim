import numpy as np
import random
from scipy.optimize import newton

class ellipticBeam:
    
    def __init__(self, _t, _c, _beta, _betaPrime=[0.,0.]):
        """ Generate a matched bunch for a fixed emittance
        Args:
        t (float) the elliptic potential strength
        c (float) the elliptic potential c
        beta (array) the beta function where the bunch is being matched
        betaPrime (array) the derivative of the beta function, defaults to zero
        """
        self.ellipticT  = -1.*_t
        self.ellipticC  = _c
        self.betaX      = _beta[0]
        self.betaY      = _beta[1]
        self.betaPrimeX = _betaPrime[0]
        self.betaPrimeY = _betaPrime[1]

    def computeHamiltonian(self, xHat, pxHat, yHat, pyHat):
        """Compute the Hamiltonian (1st invariant) for the integrable elliptic potential"""

        quadratic = 0.5 * (pxHat**2 + pyHat**2) #+ 0.5 * (xHat**2 + yHat**2)

        elliptic = 0.
        kfac = 1.
        if self.ellipticT != 0.:
            xN = xHat / self.ellipticC
            yN = yHat / self.ellipticC

            # Elliptic coordinates
            u = ( np.sqrt((xN + 1.)**2 + yN**2) +
                  np.sqrt((xN - 1.)**2 + yN**2) )/2.
            v = ( np.sqrt((xN + 1.)**2 + yN**2) -
                  np.sqrt((xN - 1.)**2 + yN**2) )/2.

            f2u = u * np.sqrt(u**2 - 1.) * np.arccosh(u)
            g2v = v * np.sqrt(1. - v**2) * (-np.pi/2 + np.arccos(v))

            kfac = self.ellipticT * self.ellipticC**2
            elliptic = (f2u + g2v) / (u**2 - v**2)

        hamiltonian = quadratic + self.computePotential(xHat, yHat)
        return hamiltonian
        
    def computePotential(self, xHat, yHat):
        quadratic = 0.5 * (xHat**2 + yHat**2)

        elliptic = 0.
        kfac = 1.
        if self.ellipticT != 0.:
            xN = xHat / self.ellipticC
            yN = yHat / self.ellipticC

            # Elliptic coordinates
            u = ( np.sqrt((xN + 1.)**2 + yN**2) +
                  np.sqrt((xN - 1.)**2 + yN**2) )/2.
            v = ( np.sqrt((xN + 1.)**2 + yN**2) -
                  np.sqrt((xN - 1.)**2 + yN**2) )/2.

            f2u = u * np.sqrt(u**2 - 1.) * np.arccosh(u)
            g2v = v * np.sqrt(1. - v**2) * (-np.pi/2 + np.arccos(v))

            kfac = self.ellipticT * self.ellipticC**2
            elliptic = (f2u + g2v) / (u**2 - v**2)

        potential = quadratic + kfac * elliptic
        return potential
        
    def whatsLeft(self, yHat):
        return self.emittance - self.computePotential(0, yHat)

    def generateBunch(self, emittance, nParticles):
        """ Generate a matched bunch with single emittance and number of particles
        Args:
        emittance (float) the value of fixed H
        nParticles(int)   the number of particles for the bunch
        
        Returns:
        bunch (list)  a list of numpy arrays of 4D phase space, (x, px, y, py)
        """
        
        # Generate some bounds on the transverse size to reduce waste in generating the bunch
        
        # Use the lemming method to find the maximum y
        y0 = np.sqrt(emittance)
        #dy = 0.01*self.ellipticC
        #while self.computePotential(0, y0) < emittance:
        #    print self.computePotential(0,y0)
        #    y0 += dy
        
        #yMax = y0
        self.emittance = emittance
        yMax = newton(self.whatsLeft, y0)        
        
        # x is harder to bound due to the peanut nature of the potential -- estimate using conventional elliptic bunch
        xMax = yMax
        
        
        # Generate particles by creating trials and finding particles with potential less than emittance, then assign the rest to momentum
        ptclsMade = 0
        phaseSpaceList = []
        while ptclsMade < nParticles:
            xTrial = 2.*(0.5 - random.random())*xMax
            yTrial = 2.*(0.5 - random.random())*yMax
            trialValue = self.computePotential(xTrial, yTrial)
            if trialValue < emittance:
                pMag = np.sqrt(2*(emittance - trialValue))
                pDir = 2*np.pi * random.random()
                pxHat = pMag * np.cos(pDir)
                pyHat = pMag * np.sin(pDir)
                xReal = xTrial * np.sqrt(self.betaX)
                yReal = yTrial * np.sqrt(self.betaY)
                pxReal = (pxHat + 0.5*self.betaPrimeX*xTrial)/np.sqrt(self.betaX)
                pyReal = (pyHat + 0.5*self.betaPrimeY*yTrial)/np.sqrt(self.betaY)
                ptclCoords = np.array([xReal, pxReal, yReal, pyReal])
                phaseSpaceList.append(ptclCoords)
                ptclsMade += 1        
        
        return phaseSpaceList
