PtrList<volScalarField> Lambda(Y.size());
forAll(Y,i)
{
word name("Lambda"+Y[i].name());
Lambda.set
(
i,new volScalarField
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
    dimensionedScalar(name,dimensionSet (1,1,-3,-1,0,0,0),0)
    )
);
}


//**********************lambda***********************//

List<dimensionedScalar> Lambda1(Y.size());
List<dimensionedScalar> Lambda2(Y.size());
List<dimensionedScalar> Lambda3(Y.size());
List<dimensionedScalar> Lambda4(Y.size());
//**********************lambda**************************//
volScalarField mixLambdaA
(
    IOobject
    (
        "mixLambdaA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("mixLambdaA",dimensionSet (1,1,-3,-1,0,0,0),0.0)
);

volScalarField Lambdai
(
    IOobject
    (
        "Lambdai",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Lambdai",dimensionSet (1,1,-3,-1,0,0,0),0.0)
);

volScalarField mixLambdaB
(
    IOobject
    (
        "mixLambdaB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("mixLambdaB",dimensionSet (-1,-1,3,1,0,0,0),0.0)
);

volScalarField mixAlpha
(
    IOobject
    (
        "mixAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mixAlpha",dimensionSet (1,-1,-1,0,0,0,0),0.0)
);
