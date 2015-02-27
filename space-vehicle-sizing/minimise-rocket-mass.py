# Speculative optimiser for liquid-fuelled suborbital space vehicle!
# The current version is kinda bad, in that it optimises for minimum mass rather
# than minimum cost. Mostly that's because I don't know how much things cost.

# TODO: learn how much things cost.

## Imports
from math import *
import numpy as np

## Handy Classes
class gas(object):
    def __init__(self,gamma,M,Cstar,rho):
        self.gamma=gamma
        self.M=M
        self.Cstar=Cstar
        self.rho=rho

class engine(object):
    def __init__(self,exhaustGas,chamberPressure,exitPressure,throatDiameter=None,seaLevelThrust=None):
        if throatDiameter==None:
            print("thrust is something!");
        else:
            print("Thrust is something else");

## Handy Functions

def TankMass(wallDensity, allowableStress, diameter, midsectionLength, pressure):
    if midsectionLength>0:
        return 2*pi*(diameter/2)**2*(diameter/2+midsectionLength)*pressure*wallDensity/allowableStress;   
    else:
        return 1.5*pressure*4/3*pi*(diameter/2)**3*wallDensity/allowableStress;
    
def airDensity(altitude):
    #scale-height model
    rhoSL=1.2;
    meanTemp=260.;
    scaleHeightConstant=29.26;
    scaleHeight=meanTemp*scaleHeightConstant;
    
    return rhoSL*exp(-altitude/scaleHeight);


## Constants
# Lots of stuff is currently approximated to zeroth order. May change in later versions.

engineTMR=400;          # This is Thrust-to-Mass ratio, hence why it's so huge.
payloadMass=40;         # Useful payload
parasiticMass=60;       # Probably the wrong term. Probably overly optimistic.

tankDensity=2700;       # Aluminium
allowableTankStress=1e8;# TODO: read about fatigue limits for Al alloys


## Constants That Probably Shouldn't Be Constants

pressurantDensity   = 1.123     # kg/m³ at 1 bar, 300K
oxidiserDensity     = 1141      # LOX
fuelDensity         = 786       # IPA

contractionRatio    = 3         # unusually this is the diameter ratio.

OFRatio             = 1.781 

exhaustGas = gas(
    gamma   = 1.156,    # at throat
    M       = 22.7321,  # molecular weight at throat
    Cstar   = 1121.5,   # at throat
    rho     = 0.0517)   # exhaust density at stagnation, at a pressure of 1 bar

