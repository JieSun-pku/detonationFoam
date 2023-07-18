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
    Foam::mixtureAverage

Description
    mixture averaged transport. Mass conservation enforced by the correction 
    velocity approach.
SourceFiles
    mixtureAverage.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureAverage_H
#define mixtureAverage_H

#include "laminarTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mixtureAverage Declaration
\*---------------------------------------------------------------------------*/

class mixtureAverage
:
    public laminarTransport
{
    // Protected data
protected:
        //-the mixture-averaged diffusion coefficients
        PtrList<volScalarField> Dmix_;
        //-correction velocity to enforce mass conservation
        mutable volVectorField Vcorr_;
        mutable volScalarField DInv;
        //- flux due to the correction velocity
        mutable surfaceScalarField phiCorr_;
        //-flux due to thermophoretic effect
        mutable surfaceScalarField phiT_;
        //--flux due to th
        mutable surfaceScalarField phiM_;
        //-thermal conductivity of species        
        PtrList<volScalarField> kappaSpecies_;
        //-temporary variables to compute mixture thermal conductivity
        mutable volScalarField kappaTemp_;
        mutable volScalarField kappaInv_;
        dictionary mixtureAverageDict_;
        //Conductivity coefficients
        scalarField b0_;
        scalarField b1_;
        scalarField b2_;
        scalarField b3_;
        //thermophoretic coefficients
        simpleMatrix<scalar> a_;
        simpleMatrix<scalar> b_;
        simpleMatrix<scalar> c_;
        simpleMatrix<scalar> d_;
        PtrList<volScalarField> Theta_;
        Switch gradX_;
        Switch thermophoresis_;
        const scalar CutOff_;
 

private:
    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixtureAverage(const mixtureAverage&);

        //- Disallow default bitwise assignment
        void operator=(const mixtureAverage&);
        //-phi, used to calculate the mixture viscosity using Wilke's formula
        tmp<volScalarField> phi(const label specieI) const;
        //-calculate the correction velocity
        void correct();
        //-access to the thermal diffusion velocity of a species
        tmp<volVectorField> VT(const label specieI) const;
        tmp<volVectorField> VT(const word& specieName) const;

public:

    //- Runtime type information
    TypeName("mixtureAverage");
    // Constructors
    mixtureAverage
    (
        const volVectorField& U,
        basicMultiComponentMixture& composition,
        basicSpecieMixture& tempcomposition,
        moleFraction& moleFraction_,
        psiReactionThermo& thermo,
        const fvMesh& mesh
    );






    //- Destructor
    ~mixtureAverage();


    // Member Functions

        // Access
        //- Return the source term for the mass fraction equations
        virtual tmp<fvScalarMatrix> Yflux(const volScalarField& Yi) const;
        
        // Update
        virtual void update();

        // Write
        virtual void write();
       



};

} // End namespace Foam

#endif

// ************************************************************************* //
