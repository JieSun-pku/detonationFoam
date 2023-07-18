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

#include "moleFraction.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(moleFraction, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleFraction::moleFraction
(
    basicMultiComponentMixture& composition,
    basicSpecieMixture& tempcomposition,
    const fvMesh& mesh
)
:
    species_(composition.species()),
    W_(species_.size()),
    X_(species_.size()),
    Y_(composition.Y()),
    sum_
    (
        IOobject
        (
            "sum",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("zero", dimless, 0.0)
    )
    
{
    word preSpecie="X.";
    forAll(X_, i)
    {
        // No need to check if field exists and can be read
        word SpecieMoleFrac = preSpecie + Y_[i].name();
        X_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    SpecieMoleFrac,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
        W_[i] = tempcomposition.W(i);

    }
    
    //- initialize the mole fractions
    sum_ *= 0.0;
    forAll(X_, i)
    {
        X_[i] = Y_[i]/W_[i];
	sum_ += X_[i];
    }
    forAll(X_, i)
    {
        X_[i] /= sum_;
    }
    
}


// ************************************************************************* //
