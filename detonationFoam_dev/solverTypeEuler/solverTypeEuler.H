while (runTime.run())
{
#include "detoCellular.H"

runTime++;
mesh.update();

#include "calculateDt.H"

Info<< "Time = " << runTime.timeName() << nl << endl;
//fluxSchemeFields->update(rho,U,e,p,c,phi,rhoPhi,rhoUPhi,rhoEPhi);
fluxSchemeFields->update(rho,rhoU,e,rPsi,c,phi,rhoPhi,rhoUPhi,rhoEPhi);

#include "rhoEqn_Euler.H"
#include "rhoUEqn_Euler.H"
#include "rhoYEqn_Euler.H"
#include "rhoEEqn_Euler.H"

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
