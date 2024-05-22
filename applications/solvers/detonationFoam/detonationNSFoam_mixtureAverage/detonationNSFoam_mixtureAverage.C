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

Application
    detonationNSFoam

Description
    Density-based compressible reactive flow solver for Euler equations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"        
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "turbulentFluidThermoModel.H"
#include "directionInterpolate.H"
#include "fvcSmooth.H"
#include "fluxScheme.H"
#include "fvOptions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    #include "readChemistryProperties.H"
    #include "createDiffFields.H"
    #include "createTDiffFields.H"
    #include "createLambdaFields.H"
    #include "createMuFields.H"

    turbulence->validate();
    Info<< "\nStarting time loop\n" << endl;

    #include "readTranData.H"
    #include "readBinaryDiff.H"
    #include "readThermalDiff.H"
    #include "readSpeciesLambda.H"
    #include "readSpeciesMu.H"

    while (runTime.run())
    {
        #include "updateTransProperties.H"
        #include "calculateDt.H"
        #include "detoCellular.H"

        mesh.update();               

        fluxSchemeFields->update(rho,U,e,p,c,phi,rhoPhi,rhoUPhi,rhoEPhi);

        #include "rhoEqn.H"
        #include "rhoUEqn.H"
        #include "rhoYEqn.H"
        #include "rhoEEqn.H"

        turbulence->correct();

        runTime.write();
    }

    Info<< "End" << endl;

    return 0;
}
