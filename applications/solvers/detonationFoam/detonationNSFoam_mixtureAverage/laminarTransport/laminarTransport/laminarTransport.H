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
    Foam::laminarTransport

Description
    base class for laminar transport models.
SourceFiles
    laminarTransport.C

\*---------------------------------------------------------------------------*/

#ifndef laminarTransport_H
#define laminarTransport_H

#include "psiReactionThermo.H"
#include "basicMultiComponentMixture.H"
#include "basicSpecieMixture.H"
#include "simpleMatrix.H"
#include "moleFraction.H"
#include "fvCFD.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class laminarTransport Declaration
\*---------------------------------------------------------------------------*/

class laminarTransport
{
 
protected:
    
        //- logarithm of T
        volScalarField logT;
        //- temporary storage for the properties
        volScalarField property;    
        //-mixture type 
        basicMultiComponentMixture& composition_;
        //-access to thermo model
        basicSpecieMixture& tempcomposition_;
        psiReactionThermo& thermo_;
        //-number of species
        const label n_;
        //- reference to species mass and mole fractions
        PtrList<volScalarField>& Y_;
        PtrList<volScalarField>& X_;
        //- viscosity of the species
        PtrList<volScalarField> muSpecies_;
        //- storage for binary diffusion coefficients
        PtrList<volScalarField> D_;
        //- mixture viscosity 
        volScalarField mu_;
        //-mixture thermal diffusivity for enthalpy (default)
        volScalarField alpha_;
        //-mixture thermal diffusivity for internal energy
        volScalarField alphaE_;
        //-mixture thermal conductivity
        volScalarField kappa_;
        //- diffusion velocities
        PtrList<volVectorField> V_;
        //- reference to velocity field
        const volVectorField& U_;
        //Viscosity Coefficients
        scalarField a0_;
        scalarField a1_;
        scalarField a2_;
        scalarField a3_;
        //diffusivity coefficients
        simpleMatrix<scalar> d0_;
        simpleMatrix<scalar> d1_;
        simpleMatrix<scalar> d2_;
        simpleMatrix<scalar> d3_;
         
        dictionary transportDict_;
        Switch viscousDissipation_;
    // Private Member Functions

        //- Disallow default bitwise copy construct
        laminarTransport(const laminarTransport&);

        //- Disallow default bitwise assignment
        void operator=(const laminarTransport&);


public:
    
     //- Runtime type information
        TypeName("laminarTransport");

    // Constructors
        //- Construct from components
        laminarTransport
        (
            const volVectorField& U,
            basicMultiComponentMixture& composition,
            basicSpecieMixture& tempcomposition,
            moleFraction& moleFraction_,
            psiReactionThermo& thermo,
            const fvMesh& mesh
        );

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            laminarTransport,
            transportModel,
            (
                const volVectorField& U,
                basicMultiComponentMixture& composition,
                basicSpecieMixture& tempcomposition,
                moleFraction& moleFraction_,
                psiReactionThermo& thermo,
                const fvMesh& mesh
            ),
            (U, composition, tempcomposition, moleFraction_, thermo, mesh)
        );        



    // Selectors

        //- Select null constructed
        static autoPtr<laminarTransport> New
        (
            const volVectorField& U,
            basicMultiComponentMixture& composition,
            basicSpecieMixture& tempcomposition,
            moleFraction& moleFraction_,
            psiReactionThermo& thermo,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~laminarTransport(){}


    // Member Functions

        // Access
        //- Constant access to the mixture dynamic viscosity
        const volScalarField& mu() const
        {
            return mu_;
        }
        //- Return the mixture kinematic viscosity
        virtual tmp<volScalarField> nu() const;
        //- Constant access to the mixture thermal diffusivity (enthalpy)
        const volScalarField& alpha() const
        {
            return alpha_;
        }
        //- Constant access to the mixture thermal diffusivity (internal 
        //  energy) 
        const volScalarField& alphaE() const
        {
            return alphaE_;
        }
        //- Return the mixture thermal conductivity
        const volScalarField& kappa() const
        {
            return kappa_;
        }
        //- Return the diffusion flux term for the mass fraction equations
        virtual tmp<fvScalarMatrix> Yflux(const volScalarField& Yi) const = 0;
        //- Net diffusive flux. Should be zero
        tmp<surfaceScalarField> sumJ() const;
        //- Return the divergence of interdiffusion sensible enthalpy flux
        tmp<volScalarField> JHs() const;
        //- Return the explicit part of heat conduction in terms of enthalpy
        tmp<volScalarField> Hconduction() const;
        //-return the explicit part of heat conduction in terms of internal energy
        tmp<volScalarField> Econduction() const;
        //-index for accessing the binary diffusivities. Since the matrix is
        // symmetric (Dij = Dji), only the upper triangle is stored.
        const label index( const label specieI, const label specieJ) const;

        // Calculate transport coefficients
        virtual void update() = 0;
        //-update the binary diffusion coefficients
        void updateBinaryDiffCoeffs();
        // read
        void read(const fvMesh& mesh);
        // Write
        virtual void write() = 0;
        //- Return the viscous stress tensor
        virtual tmp<volSymmTensorField> rhoTau() const;
        
        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divRhoTau(volVectorField& U) const;
        //- Molecular weight of the mixture [kg/kmol]
        tmp<volScalarField> W() const;
        //-the viscous dissipation
        tmp<volScalarField> viscousDissipation() const;
 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
