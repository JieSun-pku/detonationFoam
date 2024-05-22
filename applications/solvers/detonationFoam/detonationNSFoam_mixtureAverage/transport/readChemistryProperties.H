Info<< "Reading chemistry properties\n" << endl;

IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

Switch differentialDiffusion(chemistryProperties.lookupOrDefault("differentialDiffusion",true));
Switch thermalDiffusion(chemistryProperties.lookupOrDefault("thermalDiffusion",false));
dimensionedScalar Prt(chemistryProperties.lookup("Prt"));
dimensionedScalar Sct(chemistryProperties.lookup("Sct"));