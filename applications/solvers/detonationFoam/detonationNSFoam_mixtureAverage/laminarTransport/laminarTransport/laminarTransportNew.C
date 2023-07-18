/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | Code based on OpenFOAM
  \\    /   O peration     |
   \\  /    A nd           | Copyright (C) Adhiraj Dasgupta
    \\/     M anipulation  |                     
-------------------------------------------------------------------------------
 License
     This file is a derivative work of OpenFOAM.
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


\*---------------------------------------------------------------------------*/

#include "laminarTransport.H"
#include "Switch.H"

// ************************************************************************* //

Foam::autoPtr<Foam::laminarTransport> Foam::laminarTransport::New
(
    const volVectorField& U,
    basicMultiComponentMixture& composition,
    basicSpecieMixture& tempcomposition,
    moleFraction& moleFraction_,
    psiReactionThermo& thermo,
    const fvMesh& mesh    
)
{
    word transportModel;
    {
        IOdictionary transportPropertiesDict
        (
            IOobject
            (
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        transportPropertiesDict.lookup("transportModel") >> transportModel;
    }
    Info<< "Selecting transport model " << transportModel << endl;
    
    transportModelConstructorTable::iterator cstrIter =
        transportModelConstructorTablePtr_->find(transportModel);
    if (cstrIter == transportModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "laminarTransport::New"
            "("
            "    const volVectorField& U,"
            "    basicMultiComponentMixture& composition,"
            "    basicSpecieMixture& tempcomposition,"
            "    moleFraction& moleFraction_,"
            "    psiReactionThermo& thermo,"
            "    const fvMesh& mesh"
            ")"
        ) << "Unknown transport model "
          << transportModel << endl << endl
          << "Valid transport models are : " << endl
          << transportModelConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }
    return autoPtr<laminarTransport>
    (
        cstrIter()(U, composition, tempcomposition, moleFraction_, thermo, mesh)
    );
    
}
