solve(fvm::ddt(rhoU) + fvc::div(rhoUPhi));

U.ref() =
    rhoU()
   /rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

