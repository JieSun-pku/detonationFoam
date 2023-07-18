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

#include "LewisNumber.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(LewisNumber, 0);
    addToRunTimeSelectionTable
    (
        laminarTransport, 
        LewisNumber,
        transportModel
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::LewisNumber::phi
(
    const label specieI
) const
{
    tmp<volScalarField> tphi
    (
        new volScalarField
        (
            IOobject
            (
                "phiTemp",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    volScalarField& phi = tphi.ref();
//     phi *= 0.0;
    forAll(Y_, specieJ)
    {
        phi += 
        (
            1.0/sqrt(8.0)
           *pow
            (
                (
                    dimensionedScalar("one", dimless, 1.0)
                  + sqrt(muSpecies_[specieI]/muSpecies_[specieJ])
                  * pow
                    (
                        tempcomposition_.W(specieJ)
                       /tempcomposition_.W(specieI)
                       ,
                       0.25
                    )    
                ),
                2.0
            )
           *X_[specieJ] 
           /sqrt
            (
                dimensionedScalar("one", dimless, 1.0)
              + tempcomposition_.W(specieI)/tempcomposition_.W(specieJ)  
            )
        );
    }
    
    return tphi;
}

void Foam::LewisNumber::correct()
{
    Vcorr_ = 
    (
        dimensionedVector
        (
            "zero",
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            vector(0, 0, 0)
        )
    );
    forAll(Y_, i)
    {
        Vcorr_ -= Y_[i]*V_[i];
    }
    phiCorr_ = 
        linearInterpolate
        (
            thermo_.rho()
           *Vcorr_            
        )
      & kappa_.mesh().Sf();
    
    forAll(Y_, specieI)
    {
        V_[specieI] += Vcorr_;
    }
}
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LewisNumber::LewisNumber
(
    const volVectorField& U,
    basicMultiComponentMixture& composition,
    basicSpecieMixture& tempcomposition,
    moleFraction& moleFraction_,
    psiReactionThermo& thermo,
    const fvMesh& mesh
)
:
    laminarTransport(U, composition, tempcomposition, moleFraction_, thermo, mesh),
    Dmix_(n_),
    Vcorr_
    (
        IOobject
        (
            "Vcorr",
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE             
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            vector(0, 0, 0)
        )
    ),    
    phiCorr_
    (
        IOobject
        (
            "phiCorr",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0 )     
   
    ),
    kappaSpecies_(n_),
    kappaTemp_
    (
        IOobject
        (
            "kappaTemp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,1,-3,-1,0,0,0), 0.0)
    ),
    kappaInv_
    (
        IOobject
        (
            "kappaInv",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, -1, 3, 1, 0, 0, 0), 0.0)
    ),    
    LewisNumberDict_(transportDict_.subDict("LewisNumberProperties")),
    b0_(n_),
    b1_(n_),
    b2_(n_),
    b3_(n_),
    Le_(n_, 1.0)
{
   forAll(kappaSpecies_, i)
   {
       const word name = "kappa." + Y_[i].name();
       kappaSpecies_.set
       (
           i,
           new volScalarField
           (
               IOobject
               (
                   name,
                   mesh.time().timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               mesh,
               dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0) 
           )          
           
       );
   }

   forAll(Dmix_, i)
   {
       const word name = "Dmix." + Y_[i].name();
       Dmix_.set
       (
           i,
           new volScalarField
           (
               IOobject
               (
                   name,
                   mesh.time().timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               mesh,
               dimensionedScalar("zero",dimensionSet(1,-1,-1,0,0,0,0),0 )       
         
           )          
           
       );
   }

    IOdictionary conductivityPropertiesDict
    (
        IOobject
        (
            "conductivityProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false            
        )
    );
    forAll(Y_, i)
    {
        dictionary CoeffsDict = 
           conductivityPropertiesDict.subDict(Y_[i].name());
        b0_[i] = readScalar(CoeffsDict.lookup("b0"));
        b1_[i] = readScalar(CoeffsDict.lookup("b1"));
        b2_[i] = readScalar(CoeffsDict.lookup("b2"));
        b3_[i] = readScalar(CoeffsDict.lookup("b3"));
    }
    
    dictionary LeDict = LewisNumberDict_.subDict("Le");
    forAll(Y_, i)
    {
        if (LeDict.found(Y_[i].name()))
        {
            Le_[i] = readScalar(LeDict.lookup(Y_[i].name()));
        }
        else
        {
            Le_[i] = LeDict.lookupOrDefault("default", 1.0);
        }
    }
    

}    


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LewisNumber::~LewisNumber()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::LewisNumber::update()
{
    Info<< "Calculating transport coefficients" << endl;
    logT = 
    log(
           thermo_.T()
          /dimensionedScalar
           (
               "zero",
               dimensionSet(0, 0, 0, 1, 0, 0, 0),
               1.0
           )
       );
    forAll(Y_, i)
    {
        property = a0_[i] + logT*(a1_[i] + logT*(a2_[i] + a3_[i]*logT));
        muSpecies_[i] = 
        (
            exp(property)
           *dimensionedScalar
            (
                "zero",
                dimensionSet(1, -1, -1, 0, 0, 0, 0),
                1.0
            )    
        );
        
        property = b0_[i] + logT*(b1_[i] + logT*(b2_[i] + b3_[i]*logT));
        kappaSpecies_[i] = 
        (
            exp(property)
           *dimensionedScalar
            (
                "zero",
                dimensionSet(1, 1, -3, -1, 0, 0, 0),
                1.0
            )    
        );        
    }
    
    //-calculate the mixture thermal conductivity
    kappaTemp_ *= 0;
    kappaInv_ *= 0;
    forAll(X_, i)
    {
        kappaTemp_ += 
            kappaSpecies_[i]
           *(
               X_[i] + dimensionedScalar("zero", dimless, SMALL)
            );
        kappaInv_ +=
        (
           (X_[i] + dimensionedScalar("zero", dimless, SMALL))
          /kappaSpecies_[i]
        );        
    }
    kappa_ = 0.5*(kappaTemp_+ 1.0/kappaInv_);
    alpha_ = kappa_/thermo_.Cp();
    alphaE_ = kappa_/thermo_.Cv();
    
    forAll(Y_, i)
    {
        Dmix_[i] = alpha_/Le_[i];
    }
    correct();
    
    //-calculate the mixture dynamic viscosity
    mu_ *= 0;
    forAll(Y_, i)
    {
        mu_ += X_[i]*muSpecies_[i]/phi(i);
    }
    
    forAll(Y_, specieI)
    {
        V_[specieI] =
        (
            -Dmix_[specieI]*fvc::grad(Y_[specieI], "grad(Yi)")
            /(
                thermo_.rho()
               *(Y_[specieI] + dimensionedScalar("zero", dimless,SMALL)) 
             )
        );
    }    
}



void Foam::LewisNumber::write()
{
    forAll(Y_, i)
    {
        muSpecies_[i].write();
        kappaSpecies_[i].write();
        Dmix_[i].write();
        V_[i].write();
    }
    Vcorr_.write();
    rhoTau()().write();
    sumJ()().write();
    viscousDissipation()().write();
}


Foam::tmp<Foam::fvScalarMatrix> Foam::LewisNumber::Yflux
(
    const volScalarField& Yi
) const
{
    const label i = composition_.species()[Yi.name()];
    tmp<fvScalarMatrix> tYflux
    (
        - fvm::laplacian(Dmix_[i], Yi, "laplacian(Di,Yi)")
        + fvm::div(phiCorr_, Yi, "div(phi,Yi)")
    );
 
    return tYflux;
}




// ************************************************************************* //
