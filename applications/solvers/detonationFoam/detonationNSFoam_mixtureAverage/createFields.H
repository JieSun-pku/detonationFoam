Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

basicSpecieMixture& tempcomposition = thermo.composition();

basicMultiComponentMixture& Xtempcomposition = thermo.composition();

moleFraction moleFraction_(Xtempcomposition, tempcomposition, mesh);
PtrList<volScalarField>& X = moleFraction_.X();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}
volScalarField& e = thermo.he();

volScalarField maxp
(
IOobject
(
    "maxp",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
),
mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);


surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity*dimVelocity*dimArea, 0.0)
    );


surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimVelocity*dimArea, 0.0)
    );

surfaceVectorField rhoUPhi
    (
        IOobject
        (
            "rhoUPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimDensity*sqr(dimVelocity)*dimArea, Zero)
    );

surfaceScalarField rhoEPhi
    (
        IOobject
        (
            "rhoEPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity*pow3(dimVelocity)*dimArea, 0.0)
    );

autoPtr<laminarTransport> ptransport
(
  laminarTransport::New
  (
      U,
      composition,
      tempcomposition,
      moleFraction_,
      thermo,
      mesh
  )
);
laminarTransport& transport = ptransport();

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
);

volScalarField Rdot
(
    IOobject
    (
        "Rdot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("QQdot", dimEnergy/dimVolume/dimTime, 0.0)
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);


Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
    fields.add(X[i]);
    X[i].writeOpt() = IOobject::NO_WRITE;
}
fields.add(thermo.he());

autoPtr<Foam::fluxScheme> fluxSchemeFields(Foam::fluxScheme::New(mesh));
