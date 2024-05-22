solve(fvm::ddt(rhoU) + fvc::div(rhoUPhi));

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

solve
(
    fvm::ddt(rho, U) - fvc::ddt(rho, U)
  - fvm::laplacian(muEff, U)
  - fvc::div(tauMC)
);
rhoU = rho*U;
