solve(fvm::ddt(rhoU) + fvc::div(rhoUPhi));

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

solve
(
    fvm::ddt(rho, U) - fvc::ddt(rho, U)
  - fvm::laplacian(mixMu+turbulence->mut(),U)
  - fvc::div((mixMu+turbulence->mut())*dev2(Foam::T(fvc::grad(U))))
);
rhoU = rho*U;
