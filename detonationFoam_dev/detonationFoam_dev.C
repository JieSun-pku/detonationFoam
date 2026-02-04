/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor with support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "basicSprayCloud.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "psiReactionThermo.H"
#include "multivariateScheme.H"
#include "CombustionModel.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    
    doubleScalar SW_position_temp;
    doubleScalar SW_position = 0;
    doubleScalar SW_position_limit = 0;
    
    turbulence->validate();
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    scalar CoNum = 0.0;

    SW_position_limit = solverTypeProperties.lookupOrDefault<doubleScalar>("SW_position_limit",1.0);
    Info<< "SW_position_limit = "<< SW_position_limit << endl;


    word solverType(solverTypeProperties.lookupOrDefault<word>("solverType","Euler"));
    Info<< "solverType = "<< solverType << endl;    
    
    Info<< "\nStarting time loop\n" << endl;

    if(solverType=="Euler"){
        #include "solverTypeEuler.H"
    } else if(solverType=="NS_Sutherland"){
        #include "solverTypeNS_Sutherland.H"   
    } else if(solverType=="NSSpray_Sutherland"){
        #include "solverTypeNSSpray_Sutherland.H"  
    } else if(solverType=="NS_mixtureAverage"){
        #include "solverTypeNS_mixtureAverage.H"  
    } else{
        Info<< "Error solverType: "<< solverType << "!"<< nl
            << "Please set valid solverType in solverTypeProperties:" << nl
            << "Euler" << nl
            << "NS_Sutherland"<< nl
	        << "NSSpray_Sutherland"<< nl
            << "NS_mixtureAverage"<< endl; 
    }
}

// ************************************************************************* //
