surfaceScalarField sigmaDotU
(
    "sigmaDotU",
    (
        fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
      + fvc::dotInterpolate(mesh.Sf(), tauMC)
    )
  & fvc::interpolate(U)    
);

solve
(
    fvm::ddt(rhoE)
  + fvc::div(rhoEPhi)
  - fvc::div(sigmaDotU)
  ==
    reaction->Qdot()
);

e = rhoE/rho - 0.5*magSqr(U);
e.correctBoundaryConditions();
thermo.correct();
rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField()+0.5*magSqr(U.boundaryField()));


solve
(
    fvm::ddt(rho, e) - fvc::ddt(rho, e)
  + thermophysicalTransport->divq(e)
);
thermo.correct();
rhoE = rho*(e + 0.5*magSqr(U));

p.ref() = rho()/psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

Info<< "min/max(p) = "<< min(p).value() << ", " << max(p).value() << endl;
Info<< "min/max(T) = "<< min(T).value() << ", " << max(T).value() << endl;
