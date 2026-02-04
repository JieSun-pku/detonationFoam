solve(fvm::ddt(rhoU) + fvc::div(rhoUPhi));

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

/*
(
    fvm::ddt(rho, U) - fvc::ddt(rho, U)
  - fvm::laplacian(mixMu+turbulence->mut(),U)
  - fvc::div((mixMu+turbulence->mut())*dev2(Foam::T(fvc::grad(U))))
);*/

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) - fvc::ddt(rho, U)
  - fvm::laplacian(mixMu+turbulence->mut(),U)
  - fvc::div((mixMu+turbulence->mut())*dev2(Foam::T(fvc::grad(U))))
);

UEqn.relax();

UEqn.solve("U");

rhoU = rho*U;

