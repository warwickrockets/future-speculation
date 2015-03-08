# Speculative optimiser for liquid-fuelled suborbital space vehicle!
# The current version is kinda bad, in that it optimises for minimum mass rather
# than minimum cost. Mostly that's because I don't know how much things cost.

# TODO: learn how much things cost.

from rocket import *
import scipy.optimize as spop

KARMAN_LINE = 100000; #SPAAAACE!

def buildRocket(chamberPressure,exitPressure,vehicleDiameter,TMR,targetApogee):
    #Check that the dumbass optimisation algorithm has actually chosen choked flow conditions
    k=exhaustGas.gamma;
    chokedPressureRatio=(2/(k+1))**(k/(k-1))
    outerPenalties=0;

    if exitPressure<(0.4*SLPRESSURE):
        outerPenalties+=(0.4*SLPRESSURE-exitPressure)*10
        print("Separating flow is UNACCEPTABLE! "+str((0.4*SLPRESSURE-exitPressure)*10)+" YEARS DUNGEON")
        exitPressure=0.4*SLPRESSURE
        
    if chamberPressure<(SLPRESSURE/chokedPressureRatio):
        outerPenalties+=((SLPRESSURE/chokedPressureRatio)-chamberPressure)
        print("An engine isn't a vacuum chamber! "+str(((SLPRESSURE/chokedPressureRatio)-chamberPressure))+" YEARS DUNGEON")
        chamberPressure=(SLPRESSURE/chokedPressureRatio)
    
    if (exitPressure/chamberPressure > chokedPressureRatio):
        outerPenalties+=(exitPressure/chamberPressure - chokedPressureRatio)*10000
        print("Unchoked flow is UNACCEPTABLE! "+str((exitPressure/chamberPressure - chokedPressureRatio)*10000)+" YEARS DUNGEON")
        chamberPressure=exitPressure/chokedPressureRatio

    if (vehicleDiameter<0.1):
        outerPenalties+=1000*(0.1-vehicleDiameter);
        print('Such a small vehicle diameter is... UNACCEPTABLE! '+str(1000*(0.1-vehicleDiameter))+' YEARS DUNGEON');
        vehicleDiameter=0.01
        
    if (TMR<0):
        outerPenalties+=-TMR*10000
        print('You are having a bad problem and will not go to space today! Penalty: '+str(-TMR*10000)+' YEARS DUNGEON');
        TMR=-TMR
       
    print('using '+str(np.array([chamberPressure,exitPressure,vehicleDiameter,TMR])))
    
    #Initially guess that the rocket masses 400kg
    testEngine=engine(exhaustGas,chamberPressure,exitPressure,contractionRatio,seaLevelThrust=(400*TMR/engineEfficiency))
    testRocket=candidateRocket(testEngine,0.2,vehicleDiameter,0.7,300,300)

    testRocket.penalties+=outerPenalties
    
    testRocket.sizeEngineToTMR(TMR)
    
    if (targetApogee != None):
        testRocket.resizeToApogee(100,targetApogee,1)
        
    testRocket.validate();
    
    return testRocket
    
def rocketScore(theRocket):
    return theRocket.penalties+theRocket.totalLiftoffMass
    
def scoreFromFourVector(v):
    print('------')
    print('trying '+str(v))
    chamberPressure=v[0]
    exitPressure=v[1]
    vehicleDiameter=v[2]
    TMR=v[3]
    theRocket=buildRocket(chamberPressure,exitPressure,vehicleDiameter,TMR,KARMAN_LINE)
    
    return rocketScore(theRocket)
    
initialGuess=np.array([1e6,8e4,0.3,20])
spop.minimize(scoreFromFourVector,initialGuess,method="Nelder-Mead")

#stupid=buildRocket(1.00187296e+06, 7.96115084e+08, 5.57025226e+03, 9.20496522e+04, None)
#tee=np.concatenate([np.linspace(0,burnTime),np.linspace(burnTime*1.02,300,num=100)])
#y0=np.array([0,0,stupid.totalLiftoffMass])

#stupidTwo=buildRocket(7.27e4, 4.164e4,1.27189,6.20302796e+01,None)

#When resized, this guy has a terminal velocity of 7000 m/s, which means it wants more landing propellant than there is available.
#stupidThree=buildRocket(1.81580227e+05, 4.05300000e+04, 1.19586946e+00, 6.89228057e+01,100000) 
#burnTime=(stupidThree.fuelMass+stupidThree.oxMass)/stupidThree.engine.mdot
#tee=np.concatenate([np.linspace(0,burnTime),np.linspace(burnTime*1.02,300,num=100)])
#y=stupidThree.getTrajectory(tee,np.array([0.0,0.0,stupidThree.totalLiftoffMass]))