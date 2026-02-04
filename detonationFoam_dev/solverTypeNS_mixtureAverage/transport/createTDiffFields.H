//**********************thermal Diffusion****************//
List<dimensionedScalar> TDiff1(Y.size()*2);
List<dimensionedScalar> TDiff2(Y.size()*2);
List<dimensionedScalar> TDiff3(Y.size()*2);
List<dimensionedScalar> TDiff4(Y.size()*2);
//******************************************************//
volScalarField TDRatio_H
(
           IOobject
           (
              "TDRatio_H",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedScalar("TDRatio_H",dimensionSet (0,0,0,0,0,0,0),0.0)
);
volScalarField TDRatio_H2
(
           IOobject
           (
              "TDRatio_H2",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedScalar("TDRatio_H2",dimensionSet (0,0,0,0,0,0,0),0.0)
);

