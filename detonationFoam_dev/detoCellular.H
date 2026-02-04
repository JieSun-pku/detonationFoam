magGradrho = mag(fvc::grad(rho));

for (label cellI=0; cellI<mesh.C().size(); cellI++)
{
    if(maxp[cellI]<p[cellI])  maxp[cellI] = p[cellI] ;            
}

 adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
 maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
 maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", great);
 
 if (adjustTimeStep)
 {
     scalar maxDeltaTFact = maxCo/(CoNum + small);
     scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
     runTime.setDeltaT(min(deltaTFact*runTime.deltaTValue(), maxDeltaT));
     Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
 }
