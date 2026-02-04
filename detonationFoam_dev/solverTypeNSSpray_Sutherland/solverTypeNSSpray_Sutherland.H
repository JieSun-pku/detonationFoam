#include "readGravitationalAcceleration.H"
SLGThermo slgThermo(mesh, thermo);
Info<< "\nConstructing reacting cloud" << endl;
basicSprayCloud parcels
(
	 "sprayCloud",
	 rho,
	 U,
	 g,
	 slgThermo
);
#include "createRadiationModel.H"

while (runTime.run())
{
#include "detoCellular.H"
parcels.storeGlobalPositions();

runTime++;
mesh.update();

#include "calculateDt.H"

Info<< "Time = " << runTime.timeName() << nl << endl;
Info<< "Time = " << runTime.timeName() << nl << endl;
if (parcels.solution().active())
{
        parcels.evolve();
}

//fluxSchemeFields->update(rho,U,e,p,c,phi,rhoPhi,rhoUPhi,rhoEPhi);
fluxSchemeFields->update(rho,rhoU,e,rPsi,c,phi,rhoPhi,rhoUPhi,rhoEPhi);

#include "rhoEqn_NSSpray_Sutherland.H"
#include "rhoUEqn_NSSpray_Sutherland.H"
#include "rhoYEqn_NSSpray_Sutherland.H"
#include "rhoEEqn_NSSpray_Sutherland.H"

turbulence->correct();
thermophysicalTransport->correct();

runTime.write();

for (label cellI=0; cellI<mesh.C().size(); cellI++)
{
        if (p[cellI] > 101325+101325)
        {
                SW_position_temp = mesh.C()[cellI].x();
                if (SW_position_temp > SW_position)
                {
                        SW_position = SW_position_temp;
                }
        }
}
reduce(SW_position,maxOp<scalar>());
Info<< "Leading shock location = "<< SW_position << " m" << nl << endl;

if (SW_position >= SW_position_limit) exit(EXIT_SUCCESS);

Info<< "ExecutionTime = "<< runTime.elapsedCpuTime() << " s" << nl << endl;
Info<< "ClockTime = " << runTime.elapsedClockTime() << " s"<< nl << endl;
}

Info<< "End\n" << endl;

return 0;
