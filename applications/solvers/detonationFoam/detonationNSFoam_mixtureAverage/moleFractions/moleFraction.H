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

Class
    Foam::moleFraction

Description
    
    Provides a list of mole fraction fields and helper functions.

SourceFiles
    moleFraction.C

\*---------------------------------------------------------------------------*/

#ifndef moleFraction_H
#define moleFraction_H


#include "basicMultiComponentMixture.H"
#include "basicSpecieMixture.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class moleFraction Declaration
\*---------------------------------------------------------------------------*/

class moleFraction
{

protected:

    // Protected data
        //- Table of specie names
        speciesTable species_;
	//- Molecular weights
	scalarField W_;
        //- Species mole fractions
        PtrList<volScalarField> X_;
    	PtrList<volScalarField>& Y_;
	volScalarField sum_;


public:

    //- Run time type information
    TypeName("moleFraction");





    // Constructors

        //- Construct from basicMultiComponentMixture
        moleFraction
        (
           basicMultiComponentMixture&,
           basicSpecieMixture&,
	   const fvMesh&
        );


    //- Destructor
    virtual ~moleFraction()
    {}


    // Member functions

        //- Return the mole-fraction fields
        inline PtrList<volScalarField>& X();

        //- Return the const mole-fraction fields
        inline const PtrList<volScalarField>& X() const;

        //- Return the mole-fraction field for a specie given by index
        inline volScalarField& X(const label i);

        //- Return the const mole-fraction field for a specie given by index
        inline const volScalarField& X(const label i) const;

        //- Return the mole-fraction field for a specie given by name
        inline volScalarField& X(const word& specieName);

        //- Return the const mole-fraction field for a specie given by name
        inline const volScalarField& X(const word& specieName) const;
	
	//- Update the mole fractions
	inline void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "moleFractionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
