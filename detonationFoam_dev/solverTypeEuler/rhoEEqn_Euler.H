solve
(
    fvm::ddt(rhoE)
  + fvc::div(rhoEPhi)
  ==
    reaction->Qdot()
);

e = rhoE/rho - 0.5*magSqr(U);
e.correctBoundaryConditions();
thermo.correct();

//RDE
gamma=thermo.gamma();
R_gas=thermo.Cp()-thermo.Cv();
//RDE

rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField()+0.5*magSqr(U.boundaryField()));

p.ref() = rho()/psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

Info<< "min/max(p) = "<< min(p).value() << ", " << max(p).value() << endl;
Info<< "min/max(T) = "<< min(T).value() << ", " << max(T).value() << endl;
