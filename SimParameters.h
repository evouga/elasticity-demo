#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters()
    {
        timeStep = 0.0001;
        
        gravityEnabled = true;
        gravityG = -9.8;      
        lameAlpha = 100;
        lameBeta = 100;
        density = 1.0;
    }

    double timeStep;
    
    bool gravityEnabled;
    double gravityG;    
    double lameAlpha;
    double lameBeta;
    double density;
};

#endif