PtrList<volScalarField> Mu;
forAll(Y,i)
{
word name("Mu"+Y[i].name());
Mu.append
(
new volScalarField
    (
    IOobject
    (
        name,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,   
    dimensionedScalar(name,dimensionSet (1,-1,-1,0,0,0,0),0)
    )
);
}

//***********************viscosity*************************//

List<dimensionedScalar> Mu1(Y.size());
List<dimensionedScalar> Mu2(Y.size());
List<dimensionedScalar> Mu3(Y.size());
List<dimensionedScalar> Mu4(Y.size());


//********************viscocity********************//
volScalarField mixMu
(
    IOobject
    (
        "mixMu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mixMu",dimensionSet (1,-1,-1,0,0,0,0),0.0)
);

volScalarField yhi
(
    IOobject
    (
        "yhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("yhi",dimensionSet (0,0,0,0,0,0,0),0.0)
);


