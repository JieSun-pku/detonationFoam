/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    centralCourantNo

Description
    Calculates the mean and maximum wave speed based Courant Numbers.

\*---------------------------------------------------------------------------*/
bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", true);
scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
scalar maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", 1e-7);
scalar CoNum = 0.0;
dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

volScalarField rPsi("rPsi", 1.0/psi);
volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));

surfaceVectorField U_pos(interpolate(U, pos, U.name()));
surfaceVectorField U_neg(interpolate(U, neg, U.name()));

surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

surfaceScalarField cSf_pos("cSf_pos", interpolate(c, pos, T.name())*mesh.magSf());
surfaceScalarField cSf_neg("cSf_neg", interpolate(c, neg, T.name())*mesh.magSf());

surfaceScalarField ap_ref("ap_ref", max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero));
surfaceScalarField am_ref("am_ref", min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero));

surfaceScalarField a_pos("a_pos", ap_ref/(ap_ref - am_ref));
surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

surfaceScalarField amaxSf("amaxSf", max(mag(am_ref), mag(ap_ref)));
surfaceScalarField aSf_ref("aSf_ref", am_ref*a_pos);

phiv_pos *= a_pos;
phiv_neg *= a_neg;

surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf_ref);
surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf_ref);

amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

if (mesh.nInternalFaces())
{
    scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    CoNum = 0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();
}

Info<< "Max Courant Numbers = "<< CoNum << endl;

if (adjustTimeStep)
{
    scalar maxDeltaTFact = maxCo/(CoNum + small);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
 
    runTime.setDeltaT(min(deltaTFact*runTime.deltaTValue(), maxDeltaT));
 
    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

runTime++;
Info<< "Time = " << runTime.timeName() << nl << endl;

// ************************************************************************* //
