/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | Code based on OpenFOAM
  \\    /   O peration     |
   \\  /    A nd           | Copyright (C) 2017 Adhiraj Dasgupta
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

#include "mixtureAverage.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(mixtureAverage, 0);
    addToRunTimeSelectionTable
    (
        laminarTransport, 
        mixtureAverage,
        transportModel
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::mixtureAverage::phi
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

void Foam::mixtureAverage::correct()
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



Foam::tmp<Foam::volVectorField> Foam::mixtureAverage::VT
(
    const label specieI
) const
{
    const word specieName = Y_[specieI].name();
    tmp<volVectorField> tVT
    (
        new volVectorField
        (
            IOobject
            (
                "VT." + specieName,
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedVector
            (
                "zero",
                dimensionSet(0, 1, -1, 0, 0, 0, 0),
                vector(0, 0, 0)
            )
        )    
    );
    
    volVectorField& VT = tVT.ref();
    
    VT = 
    (
        -Dmix_[specieI]*Theta_[specieI]
        /(
            (X_[specieI] + dimensionedScalar("zero", dimless, SMALL))
            *thermo_.rho()
        )
        *fvc::grad(logT, "grad(T)")
    );
    
    return tVT;
    
    
}

Foam::tmp<Foam::volVectorField> Foam::mixtureAverage::VT
(
    const word& specieName
) const
{
    const label specieI = composition_.species()[specieName];
    tmp<volVectorField> tVT
    (
        new volVectorField
        (
            IOobject
            (
                "VT." + specieName,
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedVector
            (
                "zero",
                dimensionSet(0, 1, -1, 0, 0, 0, 0),
                vector(0, 0, 0)
            )
        )    
    );
    
    volVectorField& VT = tVT.ref();
    
    VT = 
    (
        -Dmix_[specieI]*Theta_[specieI]
        /(
            (X_[specieI] + dimensionedScalar("zero", dimless, SMALL))
            *thermo_.rho()
        )
        *fvc::grad(logT, "grad(T)")
    );
    
    return tVT;
    
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureAverage::mixtureAverage
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
    DInv
    (
        IOobject
        (
            "DInv",
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 1, 1, 0, 0, 0, 0), 0.0) 
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
    phiT_
    (
        IOobject
        (
            "phiT",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0 )     
   
    ),
    phiM_
    (
        IOobject
        (
            "phiM",
         mesh.time().timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0)
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, -1, 3, 1, 0, 0, 0), 0.0)
    ),    
    mixtureAverageDict_(transportDict_.subDict("mixtureAverageProperties")),
    b0_(n_),
    b1_(n_),
    b2_(n_),
    b3_(n_),
    a_(n_),
    b_(n_),
    c_(n_),
    d_(n_),
    Theta_(n_),
    gradX_(mixtureAverageDict_.lookupOrDefault("gradX", true)),
    thermophoresis_
    (
        mixtureAverageDict_.lookupOrDefault("thermophoresis", false)
    ),
    CutOff_(mixtureAverageDict_.lookupOrDefault("CutOff", 5.0))
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

   if (thermophoresis_)
   {
      Info<< "Thermophoretic effects will be considered for species with"
          << token::SPACE
          << "a molecular weight below"
          << token::SPACE
          << CutOff_
          << endl;
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
   forAll(Theta_, i)
   {
       const word name = "Theta." + Y_[i].name();
       Theta_.set
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
               dimensionedScalar("zero", dimless, 0.0)
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
   IOdictionary thermophoreticPropertiesDict
   (
       IOobject
       (
           "thermophoreticProperties",
            mesh.time().constant(),
           mesh,
           IOobject::MUST_READ,
           IOobject::NO_WRITE,
           false
       )
   );
   if (thermophoresis_)
   {
       forAll(Y_, specieI)
       {
           if (tempcomposition_.W(specieI) <= this->CutOff_)
           {
               dictionary CoeffsDict =
                   thermophoreticPropertiesDict.subDict(Y_[specieI].name());
               forAll(Y_, specieJ) 
               {
                   dictionary Coeffsij = CoeffsDict.subDict(Y_[specieJ].name());
                   a_[specieI][specieJ] = readScalar(Coeffsij.lookup("a"));
                   b_[specieI][specieJ] = readScalar(Coeffsij.lookup("b"));
                   c_[specieI][specieJ] = readScalar(Coeffsij.lookup("c"));
                   d_[specieI][specieJ] = readScalar(Coeffsij.lookup("d"));
               }
           }
       }
   }
   if(gradX_)
   {
       Info << "Diffusion velocities will be computed from mole fractions"
            << endl;
   }
   else
   {
       Info << "Diffusion velocities will be computed from mass fractions"
            << endl;
   }
}    


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixtureAverage::~mixtureAverage()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::mixtureAverage::update()
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
    updateBinaryDiffCoeffs();
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
    forAll(Y_, specieI)
    {
        DInv = 
            dimensionedScalar("zero", dimensionSet(-1, 1, 1, 0, 0, 0, 0), 0.0);
        forAll(Y_, specieJ)
        {
            if (specieJ == specieI)
            {
                continue;
            }
            const label k = index(specieI, specieJ);
            DInv +=
            (
                (X_[specieJ] + dimensionedScalar("zero", dimless, SMALL))
               /this->D_[k] 
            );
        }
        if (gradX_)
        {
            Dmix_[specieI] = 
            (
                thermo_.rho()
               /dimensionedScalar("zero",dimensionSet(1, -3, 0, 0, 0, 0, 0),1.0)
               *(dimensionedScalar("one", dimless, 1.0) - Y_[specieI])
               /DInv
            );      
        }
        else
        {
            Dmix_[specieI] = 
            (
                thermo_.rho()
               /dimensionedScalar("zero",dimensionSet(1, -3, 0, 0, 0, 0, 0),1.0)
               *(dimensionedScalar("one", dimless, 1.0) - X_[specieI])
               /DInv
            );      
        }
    }  
    if (thermophoresis_)
    {
        volScalarField Tnd =
        (
            thermo_.T()
           /dimensionedScalar
            (
               "zero",
               dimensionSet(0, 0, 0, 1, 0, 0, 0),
               1.0
            )
        );
        forAll(Y_, specieI)
        {
            Theta_[specieI] *= 0.0;
            if (tempcomposition_.W(specieI) <= this->CutOff_)
            {
                forAll(Y_, specieJ)
                {
                    Theta_[specieI] +=
                    (
                        a_[specieI][specieJ]
                      + Tnd
                       *(
                           b_[specieI][specieJ]
                         + Tnd
                          *(
                               c_[specieI][specieJ]
                             + Tnd
                              *d_[specieI][specieJ] 
                           )    
                        )   
                    )*X_[specieI]*X_[specieJ];
                }
            }
        }
    }

    forAll(Y_, specieI)
    {
        if (!gradX_)
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
        else
        {
            V_[specieI] = 
            (
                -Dmix_[specieI]*fvc::grad(X_[specieI], "grad(Xi)")
                /(
                    thermo_.rho()
                   *(X_[specieI] + dimensionedScalar("zero", dimless,SMALL))
                )
            );
        }
        V_[specieI] += VT(specieI);
    }
    correct();
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
    //-calculate the mixture dynamic viscosity
    mu_ *= 0;
    forAll(Y_, i)
    {
        mu_ += X_[i]*muSpecies_[i]/phi(i);
    }   
}



void Foam::mixtureAverage::write()
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


Foam::tmp<Foam::fvScalarMatrix> Foam::mixtureAverage::Yflux
(
    const volScalarField& Yi
) const
{
    const label i = tempcomposition_.species()[Yi.name()];
    phiT_ = dimensionedScalar("zero", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0);
    phiM_ = dimensionedScalar("zero", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0);
    if (thermophoresis_ && (tempcomposition_.W(i) <= this->CutOff_))
    {
        phiT_ = 
            linearInterpolate
            (
                thermo_.rho()
               *VT(i)
            )
            & kappa_.mesh().Sf();
    }
    if (gradX_)
    {
        phiM_ = 
            linearInterpolate
            (
               -Dmix_[i]/W()
               *fvc::grad(W()) 
            )
            & kappa_.mesh().Sf();
    }
    tmp<fvScalarMatrix> tYflux
    (
        - fvm::laplacian(Dmix_[i], Yi, "laplacian(Di,Yi)")
        + fvm::div(phiCorr_, Yi, "div(phi,Yi)")
        + fvm::div(phiT_, Yi, "div(phi,Yi)")
        + fvm::div(phiM_, Yi, "div(phi,Yi)")
    );
 
    return tYflux;
}




// ************************************************************************* //
