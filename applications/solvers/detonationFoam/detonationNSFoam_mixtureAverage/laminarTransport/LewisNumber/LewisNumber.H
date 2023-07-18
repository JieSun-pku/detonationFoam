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
    Foam::LewisNumber

Description
    Transport model based on constant Lewis number
SourceFiles
    LewisNumber.C

\*---------------------------------------------------------------------------*/

#ifndef LewisNumber_H
#define LewisNumber_H

#include "laminarTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class LewisNumber Declaration
\*---------------------------------------------------------------------------*/

class LewisNumber
:
    public laminarTransport
{
    // Protected data
protected:
        //-the mixture-averaged diffusion coefficients, computed as alpha/Lei
        PtrList<volScalarField> Dmix_;
        //-correction velocity to enforce mass conservation
        mutable volVectorField Vcorr_;
        //- flux due to the correction velocity
        mutable surfaceScalarField phiCorr_;
        //-thermal conductivity of the species        
        PtrList<volScalarField> kappaSpecies_;
        //-temporary variables to compute mixture thermal conductivity
        mutable volScalarField kappaTemp_;
        mutable volScalarField kappaInv_;
        dictionary LewisNumberDict_;
        //Conductivity coefficients
        scalarField b0_;
        scalarField b1_;
        scalarField b2_;
        scalarField b3_;
        //Lewis numbers
        scalarField Le_;
 

private:
    // Private Member Functions

        //- Disallow default bitwise copy construct
        LewisNumber(const LewisNumber&);

        //- Disallow default bitwise assignment
        void operator=(const LewisNumber&);
        //-phi, used to calculate the mixture viscosity using Wilke's formula
        tmp<volScalarField> phi(const label specieI) const;
        //-calculate the correction velocity
        void correct();
        

public:

    //- Runtime type information
    TypeName("LewisNumber");
    // Constructors
    LewisNumber
    (
        const volVectorField& U,
        basicMultiComponentMixture& composition,
        basicSpecieMixture& tempcomposition,
        moleFraction& moleFraction_,
        psiReactionThermo& thermo,
        const fvMesh& mesh
    );






    //- Destructor
    ~LewisNumber();


    // Member Functions

        // Access
        //- Return the source term for the mass fraction equations
        virtual tmp<fvScalarMatrix> Yflux(const volScalarField& Yi) const;
        
        // Update
        virtual void update();


        // Edit

        // Write
        virtual void write();
       



};

} // End namespace Foam

#endif

// ************************************************************************* //
