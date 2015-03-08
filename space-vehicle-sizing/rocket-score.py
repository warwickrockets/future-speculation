# Speculative optimiser for liquid-fuelled suborbital space vehicle!
# The current version is kinda bad, in that it optimises for minimum mass rather
# than minimum cost. Mostly that's because I don't know how much things cost.

# TODO: learn how much things cost.

from rocket import *
import scipy.optimize as spop

KARMAN_LINE = 100000; #SPAAAACE!

def buildRocket(chamberPressure,exitPressure,vehicleDiameter,TMR,targetApogee):
    
    #Initially guess that the rocket masses 400kg
    testEngine=engine(exhaustGas,chamberPressure,exitPressure,contractionRatio,seaLevelThrust=(400*TMR/engineEfficiency))
    testRocket=candidateRocket(testEngine,0.2,vehicleDiameter,0.7,300,300)
    
    
    testRocket.sizeEngineToTMR(TMR)
    
    testRocket.resizeToApogee(100,targetApogee,1)
    
    return testRocket
    
def rocketScore(theRocket):
    return theRocket.penalties-theRocket.totalLiftoffMass
    
def scoreFromFourVector(v):
    print('trying '+str(v))
    chamberPressure=v[0]
    exitPressure=v[1]
    vehicleDiameter=v[2]
    TMR=v[3]
    theRocket=buildRocket(chamberPressure,exitPressure,vehicleDiameter,TMR,KARMAN_LINE)
    
    return rocketScore(theRocket)
    
buildRocket(7.969e8,7.997e4,1.7e3,7.87e2,KARMAN_LINE)