# Tools for simulating rocket performance

## Imports
from math import *
import numpy as np
import scipy.integrate as spint
import pylab as pl

## Defs
RPRIME      = 8314;
SLPRESSURE  = 101325;
SLDENSITY   = 1.2;
G0          = 9.81;

## Handy Classes
class gas(object):
    def __init__(self,gamma,M,Cstar,T0,rho):
        self.gamma=gamma
        self.M=M
        self.Cstar=Cstar
        self.T0=T0
        self.rho=rho
    
    def densityRatio(self,machNumber):
        return ((1+(machNumber**2)*(self.gamma-1)/2)**(1/(self.gamma-1)))

class engine(object):
    def v_e(self):
        k=self.exhaustGas.gamma
        pe=self.exitPressure
        p0=self.chamberPressure
        
        return sqrt(2*self.exhaustGas.T0*(RPRIME/self.exhaustGas.M)*(k/(k-1)*(1-(pe/p0)**((k-1)/k))))

    # Specify engine either by throat diameter or sea-level thrust. Note that contractionRatio is diameter ratio!
    # See http://www.nakka-rocketry.net/th_nozz.html and http://www.nakka-rocketry.net/th_thrst.html.
    def __init__(self,exhaustGas,chamberPressure,exitPressure,contractionRatio,throatDiameter=None,seaLevelThrust=None):
        self.exhaustGas=exhaustGas
        self.chamberPressure=chamberPressure
        self.exitPressure=exitPressure

        k=exhaustGas.gamma; #typing is boring
        pe=exitPressure;
        p0=chamberPressure;

        rhoZero=self.exhaustGas.rho*chamberPressure/1e5
        rhoStar=rhoZero/self.exhaustGas.densityRatio(1);

        throatExitAreaRatio=((k+1)/2)**(1/(k-1))*(pe/p0)**(1/k)*sqrt(((k+1)/(k-1))*(1-(pe/p0)**((k-1)/k)))
        
        self.throatDiameter=0;
        
        if throatDiameter==None:
            rhoExit=rhoZero/((chamberPressure/exitPressure)**(1/k));

            A_e=seaLevelThrust/(rhoExit*self.v_e()*self.v_e()+(exitPressure-SLPRESSURE))

            self.exitDiameter=2*sqrt(A_e/pi)
            self.throatDiameter=self.exitDiameter*sqrt(throatExitAreaRatio)
        else:
            self.throatDiameter=throatDiameter
            self.exitDiameter=self.throatDiameter*sqrt(1/throatExitAreaRatio);            
        
        self.mdot=exhaustGas.Cstar*rhoStar*pi*(self.throatDiameter/2)**2
        self.chamberDiameter=self.throatDiameter*contractionRatio;
    
    def thrust(self,ambientPressure,engineEfficiency):
        return engineEfficiency*self.mdot*self.v_e()+(self.exitPressure-ambientPressure)*pi*(self.exitDiameter/2)**2;

    def v_e_effective(self,ambientPressure,engineEfficiency):
        #why would you do it this way? so ugly
        return self.thrust(ambientPressure,engineEfficiency)/self.mdot;
            
    def resizeByThrust(self,newSeaLevelThrust):
        k=self.exhaustGas.gamma; #typing is boring
        pe=self.exitPressure;
        p0=self.chamberPressure;

        rhoZero=self.exhaustGas.rho*self.chamberPressure/1e5
        rhoStar=rhoZero/self.exhaustGas.densityRatio(1);

        throatExitAreaRatio=((k+1)/2)**(1/(k-1))*(pe/p0)**(1/k)*sqrt(((k+1)/(k-1))*(1-(pe/p0)**((k-1)/k)))
        
        rhoExit=rhoZero/((self.chamberPressure/self.exitPressure)**(1/k));

        A_e=newSeaLevelThrust/(rhoExit*self.v_e()*self.v_e()+(self.exitPressure-SLPRESSURE))

        self.exitDiameter=2*sqrt(A_e/pi)
        self.throatDiameter=self.exitDiameter*sqrt(throatExitAreaRatio)        
 
        self.mdot=exhaustGas.Cstar*rhoStar*pi*(self.throatDiameter/2)**2
        self.chamberDiameter=self.throatDiameter*contractionRatio;
            

## Handy Functions

def TankMass(wallDensity, allowableStress, diameter, midsectionLength, pressure):
    if midsectionLength>0:
        return 2*pi*(diameter/2)**2*(diameter/2+midsectionLength)*pressure*wallDensity/allowableStress;   
    else:
        return 1.5*pressure*4/3*pi*(diameter/2)**3*wallDensity/allowableStress;

def airDensity(altitude):
    #scale-height model.
    #TODO: use better density model
    rhoSL=1.2;
    meanTemp=260.;
    scaleHeightConstant=29.26;
    scaleHeight=meanTemp*scaleHeightConstant;
    
    return rhoSL*exp(-altitude/scaleHeight);

def airPressure(altitude):
    #yuck.
    #TODO: use better pressure model.
    rhoSL=1.2;
    return SLPRESSURE*airDensity(altitude)/rhoSL;

## Constants
# Lots of stuff is currently approximated to zeroth order. May change in later versions.

engineTMR=400;          # This is Thrust-to-Mass ratio, hence why it's so huge.
payloadMass=40;         # Useful payload
parasiticMass=60;       # Probably the wrong term. Probably overly optimistic.

tankDensity=2700;       # Aluminium
allowableTankStress=1e8;# TODO: read about fatigue limits for Al alloys

engineEfficiency=0.9;   # Optimistic!

## Constants That Probably Shouldn't Be Constants

pressurantDensity   = 0.1604    # kg/m³ at 1 bar, 300K. (Helium)
oxidiserDensity     = 1141      # LOX
fuelDensity         = 786       # IPA

contractionRatio=2.5;   # Ratio of diameters. Likely to want to be higher for pintle engines.

OFRatio             = 1.781 

exhaustGas = gas(
    gamma   = 1.156,    # at throat
    M       = 22.7321,  # molecular weight at throat
    Cstar   = 1121.5,   # at throat
    T0      = 3207.3,   # at stagnation
    rho     = 0.08323)  # exhaust density at stagnation, at a pressure of 1 bar

## The Vehicle

class candidateRocket(object):
    def __init__(self,candidateEngine,pressureDropRatio,vehicleDiameter,C_D,pressurantStorageTemp,propellantMass):
        self.penalties=0;
        
        if(candidateEngine.exitPressure<0.4*SLPRESSURE):
            engineToUse=engine(candidateEngine.exhaustGas,candidateEngine.chamberPressure,0.4*SLPRESSURE,(candidateEngine.chamberDiameter/candidateEngine.throatDiameter),throatDiameter=candidateEngine.throatDiameter);
            self.engine=engineToUse;
            self.penalties+=0.1*(0.4*SLPRESSURE-candidateEngine.exitPressure)
        else:
            self.engine=candidateEngine
            
        self.tankPressure=self.engine.chamberPressure*(1+pressureDropRatio)
        
        self.fuelMass=propellantMass/(1+OFRatio)
        self.oxMass=self.fuelMass*OFRatio
        
        self.fuelTankVolume=self.fuelMass/fuelDensity
        self.oxTankVolume=self.oxMass/oxidiserDensity
        
        # Assuming an ideal gas, tank mass fraction is independent of pressure, so just assume it's this size.
        self.storedPressurantVolume=(self.fuelTankVolume+self.oxTankVolume)*(pressurantStorageTemp/300);
        
        self.pressurantMass=(self.fuelTankVolume+self.oxTankVolume)*pressurantDensity*(self.tankPressure/1e5);
        
        self.vehicleDiameter=vehicleDiameter;
        self.C_D=C_D
        
        self.calculateMasses();
        
        self.validate();
        
    def calculateMasses(self):
        #First calculate the volume of a spherical tank with the same diameter as the rocket. 
        #If anything needs less space than that, put it in a smaller spherical tank. Otherwise,
        #calculate the mass of a capsule-shaped tank.
        sphericalTankVolume=4/3*pi*(self.vehicleDiameter/2)**3
        
        if self.fuelTankVolume<sphericalTankVolume:
            self.fuelTankCylLength=0
            self.oxTankCylLength=0
            
            fuelTankDiameter=2*(self.fuelTankVolume*3/(4*pi))**(1/3)
            self.fuelTankMass=TankMass(tankDensity, allowableTankStress, fuelTankDiameter, 0, self.tankPressure);
        else:
            volumeInCylindricalBit=self.fuelTankVolume-sphericalTankVolume
            self.fuelTankCylLength=volumeInCylindricalBit/(pi*(self.vehicleDiameter/2)**2)
            self.fuelTankMass=TankMass(tankDensity, allowableTankStress, self.vehicleDiameter, self.fuelTankCylLength, self.tankPressure);
            
        if self.oxTankVolume<sphericalTankVolume:
            oxTankDiameter=2*(self.oxTankVolume*3/(4*pi))**(1/3)
            self.oxTankMass=TankMass(tankDensity, allowableTankStress, oxTankDiameter, 0, self.tankPressure);
        else:
            volumeInCylindricalBit=self.oxTankVolume-sphericalTankVolume
            self.oxTankCylLength=volumeInCylindricalBit/(pi*(self.vehicleDiameter/2)**2)
            self.oxTankMass=TankMass(tankDensity, allowableTankStress, self.vehicleDiameter, self.oxTankCylLength, self.tankPressure);        
        
        #Recall that 'tankage efficiency' is independent of pressure, given constant temperature.
        self.pressurantTankMass=1.5*1e5*self.storedPressurantVolume*tankDensity/allowableTankStress;
        
        self.totalTankageMass=self.pressurantTankMass+self.oxTankMass+self.fuelTankMass;
        self.engineMass=self.engine.thrust(SLPRESSURE,engineEfficiency)/engineTMR;
        self.otherDryMass=payloadMass+parasiticMass;
        
        self.finalBurnoutMass=self.totalTankageMass+self.engineMass+self.otherDryMass+self.pressurantMass;
        self.totalLiftoffMass=self.finalBurnoutMass+self.fuelMass+self.oxMass;   
        
        self.landingPropellantMass=self.finalBurnoutMass*(exp((self.emptyTerminalVelocity())/(self.engine.v_e_effective(SLPRESSURE, engineEfficiency)))-1)
        
        if self.landingPropellantMass>(0.95*(self.fuelMass+self.oxMass)):
            self.penalties+=self.landingPropellantMass-(0.95*(self.fuelMass+self.oxMass))
            print("That much landing propellant is UNACCEPTABLE! "+str(self.landingPropellantMass-(0.95*(self.fuelMass+self.oxMass)))+" YEARS DUNGEON")
            self.landingPropellantMass=0.95*(self.fuelMass+self.oxMass);
        
        self.MECOMass=self.finalBurnoutMass+self.landingPropellantMass;
        
        
    def emptyTerminalVelocity(self):
        return sqrt((2*G0*self.finalBurnoutMass)/(SLDENSITY*pi*self.C_D*(self.vehicleDiameter/2)**2));
             
    def validate(self,noPenalties=None):
        engineEnvelope=max(self.engine.exitDiameter,self.engine.chamberDiameter);
        
        if(engineEnvelope>self.vehicleDiameter):
            if noPenalties!=True:
                self.penalties+=100*((engineEnvelope-self.vehicleDiameter)/self.vehicleDiameter);
                print('This engine size is... UNACCEPTABLE! '+str(self.penalties)+' YEARS DUNGEON');
            
            self.vehicleDiameter=engineEnvelope;
            self.calculateMasses();
            

    def ydot(self, y, t):
        x=y[0]
        v=y[1]
        m=y[2]
        
        xdot=v
        
        startMECOTime=(self.totalLiftoffMass-self.MECOMass-1)/self.engine.mdot
        
        MECODuration=(((self.totalLiftoffMass-self.MECOMass)/self.engine.mdot)-startMECOTime)*2
        
        throttleLevel=np.clip((((startMECOTime+MECODuration)-t)/MECODuration),0,1)
        
        mdot=-self.engine.mdot*throttleLevel;
        T=self.engine.thrust(airPressure(x),engineEfficiency)*throttleLevel;
    
        vdot=(T-0.5*airDensity(x)*self.C_D*(pi*(0.5*self.vehicleDiameter)**2)*v**2*np.sign(v))/m-G0;

        if(x<=0):
            xdot=max(xdot,0)
            vdot=max(vdot,0)

        #Hideously hacky :(
        #TODO: some sort of event detection for spint.ode???
        #if (T==0) and (v<=0):
        #    xdot=0
        #    vdot=0
        #    mdot=0
        
        return np.array([xdot,vdot,mdot])
        
    def getTrajectory(self, t, y0):
        traj=spint.odeint(self.ydot, y0, t, mxstep=2000, atol=1e-6);
        return traj;
        
    def getApogee(self):
        burnTime=(self.fuelMass+self.oxMass)/self.engine.mdot
        
        tee=np.concatenate([np.linspace(0,burnTime*1.1,num=50),np.linspace(burnTime*1.12,300,num=20)])
        y=self.getTrajectory(tee,np.array([0.0,0.0,self.totalLiftoffMass]))
        return (np.nanmax(y,axis=0))[0]
     
    def getBoosterLength(self):
         return self.fuelTankCylLength+self.oxTankCylLength+2*(self.vehicleDiameter)
        
    def sizeEngineToTMR(self,TMR):
        self.calculateMasses()
        
        if(TMR>(engineTMR/2)):
            #And even this is ridiculous: why would you let half of your vehicle be engine? Weirdo.
            self.sizeEngineToTMR(engineTMR/2)
            self.penalties+=(TMR-engineTMR/2)*50
        else:
            diff=99999
            
            while diff>0.001:
                prevMass=self.totalLiftoffMass;
                self.engine.resizeByThrust((self.totalLiftoffMass*TMR)/engineEfficiency)
                self.calculateMasses()
                diff=abs(self.totalLiftoffMass-prevMass);
        
        
    def resizeToApogee(self,maxIterations,targetApogee,epsilon):
        currentLTMR=self.engine.thrust(SLPRESSURE,engineEfficiency)/self.totalLiftoffMass;
        i=0;
        while (i<maxIterations):
            currentApogee=self.getApogee();
        
            # Rocket is already stupidly huge, or can't get off the ground
            if (self.totalLiftoffMass>2000000 or (self.engine.thrust(SLPRESSURE,engineEfficiency)/self.MECOMass)<9.81):
                self.penalties+=abs(currentApogee-targetApogee)
                return -1
        
            if (abs(currentApogee-targetApogee)<epsilon):
                return 0
            else:
                scaleFactor = (targetApogee/currentApogee)**0.5

                self.fuelMass*=scaleFactor
                self.oxMass*=scaleFactor
        
                self.fuelTankVolume*=scaleFactor
                self.oxTankVolume*=scaleFactor
        
                # Assuming an ideal gas, tank mass fraction is independent of pressure, so just assume it's this size.
                self.storedPressurantVolume*=scaleFactor
        
                self.pressurantMass*=scaleFactor
                
                self.calculateMasses()
                self.sizeEngineToTMR(currentLTMR)
                self.validate(noPenalties=True)
                
            i+=1
            
        return -1
        
#testEngine=engine(exhaustGas,5e5,8e4,3,seaLevelThrust=15000);
#testr=candidateRocket(testEngine,0.2,0.3,0.7,300,300);
#print(testr.getApogee());