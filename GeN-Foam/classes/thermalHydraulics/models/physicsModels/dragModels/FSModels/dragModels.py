### IMPORTS

from matplotlib import pyplot as plt
import math

### CLASSES

# Base abstract class for a dragModel

class dragModel :

    def __init__(self, label) :

        self.label = label

    def fd(self, Re) :
        pass

    def fds(self, Res) :
        fds = []
        for Re in Res :
            fds.append(self.fd(Re))
        return fds

    def plot(self, Res) :
        plt.loglog(Res, self.fds(Res), label=self.label)

# Derived class for each dragModel, add your own if you want! You only need
# to implement the fd function and constructor

###############################################################################

class DarcyReynoldsPower(dragModel) :

    def __init__(self, coeff, exp, label) :

        super().__init__(label)
        self.coeff = coeff
        self.exp = exp

    def fd(self, Re) :
        return self.coeff*(Re**self.exp)

###############################################################################

class BaxiDalleDonne(dragModel) :

    def __init__(self, Dp, Dw, H, label) :

        super().__init__(label)
        '''
        Dp is pin diameter
        Dw is wire diameter
        H is wire lead len
        '''
        Pt = Dp+1.0444*Dw
        self.A = (80.0/math.sqrt(H*100))*((Pt/Dp)**1.5)
        self.B = 1.034/((Pt/Dp)**0.124)
        self.C = 29.7*((Pt/Dp)**6.9)/((H/(Dp+Dw))**2.239)

    def fd(self, Re) :

        fl = self.A/Re
        ft = 0.316*(self.B+self.C*(Re**0.086))/(Re**0.25)
        fd = fl
        if Re > 400 :
            if Re < 5000 :
                psi = min(max((Re-400.0)/4600.0, 0.0), 1.0)
                fd = math.sqrt(psi)*ft+math.sqrt(1.0-psi)*fl
            else :
                fd = ft
        return fd

###############################################################################

class Rehme(dragModel) :

    def __init__(self, Dp, Dw, H, Np, wetWrapPer, label) :

        super().__init__(label)
        '''
        Dp is pin diameter
        Dw is wire diameter
        H is wire lead len
        Np is number of pins
        wetWrapPer is wetWrapperPerimeter
        '''
        wetPinPer = Np*3.1415*(Dp+Dw)
        Pt = Dp+1.0444*Dw
        B = math.sqrt(Pt/Dp) + (7.6*(Dp+Dw)*((Pt/Dp)**2)/H)**2.16

        self.A = wetPinPer/(wetPinPer+wetWrapPer)
        self.B1 = 64.0*math.sqrt(B)
        self.B2 = 0.0816*(B**0.9335)

    def fd(self, Re) :

        return (self.A*(self.B1/Re + self.B2/(Re**0.133)))

###############################################################################

# I still need to finish implementing all the other ones

class Churchill(dragModel) :

    def __init__(self, label) :

        super().__init__(label)

        pass
        

    def fd(self, Re) :

        pass

###############################################################################

### FUNCTIONS

### USER-SELECTABLE

# Reynolds (i.e. x axis) binning. Plots are all in log-log
Re0 = 10
Re1 = 200000
N = 1000

# Models to be plotted
dragModels = []
dragModels.append(Rehme(0.006, 0.0016, 0.047, 27, 0.17524, "Rehme1"))
dragModels.append(Rehme(0.006, 0.0016, 0.047, 10, 0.17524, "Rehme2"))

### MAIN

# Setup equally spaced Re axis in log space
Re0 = max(Re0, 1e-2)
c = (((Re1-Re0)/Re0)+1)**(1.0/N)
Res = [Re0*(c**i) for i in range(N+1)]

for dm in dragModels :
    dm.plot(Res)

plt.xlabel("Re (-)")
plt.ylabel("fd (-)")
plt.legend()
plt.grid(True)
plt.show()