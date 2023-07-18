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

Application
    fitTransport

Description
    Reads the transport database and fit polynomials to the values of viscosity, 
    thermal conductivity, and binary diffusivities. The binary diffusivities are 
    fitted at 1 bar pressure, so the obtained diffusivity needs to be divided 
    by the pressure in bar.
    
    The fitting is done using the GNU Scientific Library, the latest version of
    which can be obtained here:
        <https://www.gnu.org/software/gsl/>

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "argList.H"
#include "Time.H"
#include "chemkinReader.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "labelField.H"
#include "simpleMatrix.H"
#include "molecularTransport.H"
#include "Polynomial.H"
#include "makeGraph.H"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_linalg.h>

using namespace Foam;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    const scalar RR = Foam::constant::physicoChemical::R.value()*1000;
    const word& gFormat = runTime.graphFormat(); 
    IFstream controlFile(runTime.constant()/"transportProperties");
    dictionary control(controlFile);
    const scalar Tmin(control.lookupOrDefault("Tmin", 300));
    const scalar Tmax(control.lookupOrDefault("Tmax", 3000));
    const scalar dT = 100.0;
    const label step = (Tmax-Tmin)/dT + 1;
    gsl_matrix* X;
    gsl_matrix* cov;
    gsl_vector* y;
    gsl_vector* z;
    gsl_vector* z0;
    gsl_vector* z1;
    gsl_vector* a;
    gsl_vector* b;
    gsl_vector* d;
    gsl_vector* e;
    gsl_vector* w;
    double chisqVisc, chisqCond, chisqDiff, chisqTherm;
    double rcond;
    X = gsl_matrix_alloc(step, 4);
    cov = gsl_matrix_alloc(4, 4);
    
    y = gsl_vector_alloc(step);
    z = gsl_vector_alloc(step);
    z0 = gsl_vector_alloc(step);
    z1 = gsl_vector_alloc(step);
    
    a = gsl_vector_alloc(4);
    b = gsl_vector_alloc(4);
    d = gsl_vector_alloc(4);
    e = gsl_vector_alloc(4);
    w = gsl_vector_alloc(step);
    
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(step, 4);

    IOdictionary viscosityPropertiesDict
    (
        IOobject
        (
            "viscosityProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );
    IOdictionary conductivityPropertiesDict
    (
        IOobject
        (
            "conductivityProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );
    IOdictionary diffusivityPropertiesDict
    (
        IOobject
        (
            "diffusivityProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );
    IOdictionary thermophoreticPropertiesDict
    (
        IOobject
        (
            "thermophoreticProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );
    
    fileName kineticsFile
    (
        fileName(control.lookup("CHEMKINFile")).expand()
    );
    fileName thermoFile
    (
        fileName(control.lookup("CHEMKINThermoFile")).expand()
    );    
    fileName transportFile
    (
        fileName(control.lookup("CHEMKINTransportFile")).expand() 
    );

    fileName transportFileNew
    (
        fileName(control.lookup("CHEMKINTransportFileNew")).expand() 
    );

    speciesTable species;
    chemkinReader cr(species, kineticsFile, transportFile, thermoFile, false);
    IFstream transport(transportFileNew);
    
    const label n = species.size();
    labelField molecule(n);
    scalarField LJP(n);
    scalarField SigmaLJ(n);
    scalarField Dipole(n);
    scalarField Polarizability(n);
    scalarField Zrot298(n);
    simpleMatrix<scalar> SigmaLJ_ij(n);
    simpleMatrix<scalar> LJP_ij(n);
    const scalar DipoleSI = 3.335641e-30;
    scalar Zeta;
    molecularTransport transportModel;
    scalarField T(step);
    PtrList<scalarField> mu(n);
    PtrList<scalarField> kappa(n);
    PtrList<scalarField> DTemp(n);
    PtrList<scalarField> DTTemp(n);
    
    forAll(T, k)
    {
        T[k] = Tmin + k*dT;
    }
    forAll(species, i)
    {
        mu.set(i, new scalarField(step));
        kappa.set(i, new scalarField(step));
        DTemp.set(i, new scalarField(step));
        DTTemp.set(i, new scalarField(step));
    }
    label count = 0;
    
    string line;
    Info<< "Reading transport database" << endl;
    
    
    while(transport.getLine(line))
    {
        char c = line.c_str()[0];
        
        if ((c != '!') && (c != ' ') && (c != '"'))
        {
            word name;
            IStringStream ins(line);         
            ins >> name;
            
            if (species.contains(name))
            {
                count++;
                label i = species[name];
                ins >> molecule[i];
                ins >> LJP[i];
                ins >> SigmaLJ[i];
                ins >> Dipole[i];
                ins >> Polarizability[i];
                ins >> Zrot298[i];
                              
                SigmaLJ[i] = SigmaLJ[i]*1e-10;
                Polarizability[i] = Polarizability[i]*1e-30;
                Dipole[i] = Dipole[i]*DipoleSI;
            }
        }
    }
    
    scalar DpM1, DpM2;
    scalar alpha_n, mu_p;
    label indexP, indexN;
    
    forAll(species, i)
    {
        forAll(species, j)
        {
            DpM1 = Dipole[i]/DipoleSI;
            DpM2 = Dipole[j]/DipoleSI;
            
            if 
            (
                ((DpM1 >= 1e-4) && (DpM2 >= 1e-4))
             || ((DpM1 < 1e-4) && (DpM2 < 1e-4)) 
            )
            {
                LJP_ij[i][j]= Foam::sqrt(LJP[i]*LJP[j]);
                SigmaLJ_ij[i][j] = 0.5*(SigmaLJ[i] + SigmaLJ[j]);
            }
            
            else
            {
                if ((DpM1 >= 1e-4) && (DpM2 < 1e-4))
                {
                    indexP = i;
                    indexN = j;
                }
                if ((DpM1 < 1e-4) && (DpM2 >= 1e-4))
                {
                    indexP = j;
                    indexN = i;
                }
                
                alpha_n = 
                    Polarizability[indexN]
                   /Foam::pow(SigmaLJ[indexN], 3.0);
                mu_p =
                    Dipole[indexP]
                   /Foam::sqrt
                    (
                        LJP[indexP]
                       *Foam::constant::physicoChemical::k.value()
                       *Foam::pow(SigmaLJ[indexP], 3.0)
                    );
                
                Zeta =
                    1.0 
                  + 0.25
                   *alpha_n*mu_p
                   *Foam::sqrt(LJP[indexP]/LJP[indexN]);
                SigmaLJ_ij[i][j] =
                    0.5
                   *(SigmaLJ[i] + SigmaLJ[j])
                   /Foam::pow(Zeta, (1.0/6.0));
                LJP_ij[i][j] = Zeta*Zeta*Foam::sqrt(LJP[i]*LJP[j]);   
            }
        }
    }
    
    Info<< "Read " << count << "  species" << endl;

    Info<< "The fits will be computed between the temperatures of"
        << token::SPACE
        << Tmin
        << "K"
        << token::SPACE
        << "and"
        << token::SPACE
        << Tmax
        << "K"
        << endl;
    OStringStream  noteFit;
    noteFit << "The fits were computed between the temperatures of"
             << token::SPACE
             << Tmin
             << "K"
             << token::SPACE
             << "and"
             << token::SPACE
             << Tmax
             <<"K";    
    viscosityPropertiesDict.note() = noteFit.str();
    conductivityPropertiesDict.note() = noteFit.str();
    diffusivityPropertiesDict.note() = noteFit.str();
    thermophoreticPropertiesDict.note() = noteFit.str();
    
    IOdictionary transportData
    (
        IOobject
        (
            "transportData",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    forAll(species, i)
    {
        dictionary speciesDict;        
        speciesDict.add("moleculeType", molecule[i]);
        speciesDict.add("LennardJonesParam", LJP[i]);
        speciesDict.add("CollisionDiam", SigmaLJ[i]);
        speciesDict.add("dipoleMoment", Dipole[i]);
        speciesDict.add("polarizability", Polarizability[i]);
        speciesDict.add("RotationalRelaxation", Zrot298[i]);
        
        transportData.add(species[i], speciesDict);
          
    }
    
    
    forAll(species, i)
    {
        HashPtrTable<gasHThermoPhysics>::const_iterator specieThermoIter = 
            cr.speciesThermo().find(species[i]);
            

        for (label j = 0; j < step; j++)
        {
            scalar T0 = T[j];
            const scalar logT = Foam::log(T0);
            scalar p0 = 1e5;
            const scalar omega22 =
                transportModel.CollisionIntegral22(T0/LJP[i]);                
            const scalar m = 
                (*specieThermoIter)->W()
               /(Foam::constant::physicoChemical::NA.value()*1000.0);  
            const scalar visc =
                transportModel.viscosity(m, T0, SigmaLJ[i], omega22);
                
            gsl_matrix_set(X, j, 0, 1.0);
            gsl_matrix_set(X, j, 1, logT);
            gsl_matrix_set(X, j, 2, logT*logT);
            gsl_matrix_set(X, j, 3, logT*logT*logT);
            gsl_vector_set(y, j, Foam::log(visc));
            gsl_vector_set(w, j, 1.0);
            
            const scalar omega11 =
                transportModel.CollisionIntegral11(T0/LJP[i]);
            const scalar selfdiffusion =
                transportModel.diffusivity(m, T0, SigmaLJ[i], omega11);
                
            const scalar Cp = (*specieThermoIter)->cp(p0, T0); 
            const scalar Cv = (*specieThermoIter)->cv(p0, T0);
            const scalar CvTrans = 1.5*RR;
            scalar CvRot, CvVib;
            if (molecule[i] == 1)
            {
                CvRot = RR;
                CvVib = Cv - 2.5*RR;                
            }
            else if (molecule[i] == 2)
            {
                CvRot = 1.5*RR;
                CvVib = Cv -3.0*RR;
            }
            else if (molecule[i] == 0)
            {
                CvRot = 0;
                CvVib = 0;
            }
            if (CvVib < 0.0)
            {
                CvVib = 0.0;
            }
            const scalar F = transportModel.F(LJP[i], T0);
            const scalar Zrot = Zrot298[i]*F;
            const scalar rho = (*specieThermoIter)->rho(p0, T0);
            const scalar fVib = (6.0/5.0)*(omega22/omega11);
            const scalar A = 2.5 - fVib;
            const scalar B =
                Zrot
              + (
                    2.0/Foam::constant::mathematical::pi
                   *(
                        5.0/3.0*(CvRot/RR)
                      + fVib 
                    )
                );
            const scalar fTrans = 
                2.5
               *(
                    1.0
                  - 2.0/Foam::constant::mathematical::pi
                   *(CvRot*A)
                   /(CvTrans*B)
                );
            const scalar fRot = 
                fVib
               *(
                    1.0
                  + 2.0/Foam::constant::mathematical::pi
                   *(A/B)
                );
            const scalar conductivity =
                visc/(*specieThermoIter)->W()
               *(
                    fTrans*CvTrans + fRot*CvRot + fVib*CvVib
                );
            gsl_vector_set(z, j, Foam::log(conductivity));   
        }
            
        gsl_multifit_linear(X, y, a, cov, &chisqVisc, work); 
        
        dictionary viscosityCoeffsDict;
        viscosityCoeffsDict.add("a0", gsl_vector_get(a, 0));
        viscosityCoeffsDict.add("a1", gsl_vector_get(a, 1));
        viscosityCoeffsDict.add("a2", gsl_vector_get(a, 2));
        viscosityCoeffsDict.add("a3", gsl_vector_get(a, 3));
        viscosityCoeffsDict.add("chiSquare", chisqVisc);
        viscosityPropertiesDict.add(species[i], viscosityCoeffsDict);
        
        gsl_multifit_linear(X, z, b, cov, &chisqCond, work);
        dictionary conductivityCoeffsDict;
        conductivityCoeffsDict.add("b0", gsl_vector_get(b, 0));
        conductivityCoeffsDict.add("b1", gsl_vector_get(b, 1));
        conductivityCoeffsDict.add("b2", gsl_vector_get(b, 2));
        conductivityCoeffsDict.add("b3", gsl_vector_get(b, 3));
        conductivityCoeffsDict.add("chiSquare", chisqCond);
        conductivityPropertiesDict.add(species[i], conductivityCoeffsDict);
       
        scalarField muCoeffs(4);
        scalarField kappaCoeffs(4);
        
        forAll(muCoeffs, s)
        {
            muCoeffs[s] = gsl_vector_get(a, s);
            kappaCoeffs[s] = gsl_vector_get(b, s);
        }
        Polynomial<4> muPoly(muCoeffs);
        Polynomial<4> kappaPoly(kappaCoeffs);
        
        for (label j = 0; j < step; j++)
        {
            scalar T0 = T[j];
            scalar logT = Foam::log(T0);
            mu[i][j] = Foam::exp(muPoly.value(logT));
            kappa[i][j] = Foam::exp(kappaPoly.value(logT));
        }
    }
    forAll(species, i)
    {
        dictionary diffusivityI;
        fileName DPath
        (
            runTime.rootPath()/runTime.caseName()/"transport"/"diffusivity"
        );
        mkDir(DPath);
        forAll(species, j)
        {
            HashPtrTable<gasHThermoPhysics>::const_iterator specieThermoIterI =
                cr.speciesThermo().find(species[i]);
            HashPtrTable<gasHThermoPhysics>::const_iterator specieThermoIterJ =
                cr.speciesThermo().find(species[j]);    
            for (label k = 0; k < step; k++)
            {
                scalar T0 = T[k];
                const scalar logT = Foam::log(T0);
                const scalar omega11 = 
                    transportModel.CollisionIntegral11(T0/LJP_ij[i][j]);
                const scalar mI = 
                    (*specieThermoIterI)->W()
                   /(Foam::constant::physicoChemical::NA.value()*1000.0);
                const scalar mJ = 
                    (*specieThermoIterJ)->W()
                   /(Foam::constant::physicoChemical::NA.value()*1000.0);
                const scalar m = mI*mJ/(mI + mJ);
                const scalar diffusivity = 
                    transportModel.diffusivity
                    (
                        m, T0, SigmaLJ_ij[i][j], omega11
                    );
                gsl_vector_set(z0, k, Foam::log(diffusivity));    
            }
            gsl_multifit_linear(X, z0, d, cov, &chisqDiff, work);
            dictionary diffusivityCoeffsDict;
            diffusivityCoeffsDict.add("d0", gsl_vector_get(d, 0));
            diffusivityCoeffsDict.add("d1", gsl_vector_get(d, 1));
            diffusivityCoeffsDict.add("d2", gsl_vector_get(d, 2));
            diffusivityCoeffsDict.add("d3", gsl_vector_get(d, 3));
            diffusivityCoeffsDict.add("chiSquare", chisqDiff);
                        
            diffusivityI.add(species[j], diffusivityCoeffsDict);
            
            scalarField DCoeffs(4);
            forAll(DCoeffs, s)
            {
                DCoeffs[s] = gsl_vector_get(d, s);
            }
            Polynomial<4> DPoly(DCoeffs);
            word DName = species[i] + "." + species[j];
            for (label k = 0; k < step; k++)
            {
                DTemp[j][k] = Foam::exp(DPoly.value(Foam::log(T[k])));
            }
            
            makeGraph(T, DTemp[j], DName, DPath, gFormat);            
        }
        diffusivityPropertiesDict.add(species[i], diffusivityI);
    }
    forAll(species, i)
    {
        dictionary thermophoreticI;
        fileName DTPath(runTime.rootPath()/runTime.caseName()/"transport"/"DT");
        mkDir(DTPath);
        forAll(species, j)
        {
            HashPtrTable<gasHThermoPhysics>::const_iterator specieThermoIterI =
                cr.speciesThermo().find(species[i]);
            HashPtrTable<gasHThermoPhysics>::const_iterator specieThermoIterJ =
                cr.speciesThermo().find(species[j]);
            const scalar Wi = (*specieThermoIterI)->W();
            const scalar Wj = (*specieThermoIterJ)->W();
            for (label k = 0; k < step; k++)
            {
                scalar T0 = T[k];
                const scalar omega11 = 
                    transportModel.CollisionIntegral11(T0/LJP_ij[i][j]);
                const scalar omega12 = 
                    transportModel.CollisionIntegral12(T0/LJP_ij[i][j]);
                const scalar omega13 = 
                    transportModel.CollisionIntegral13(T0/LJP_ij[i][j]);
                const scalar omega22 = 
                    transportModel.CollisionIntegral22(T0/LJP_ij[i][j]);
                    
//                 const scalar aStar = 0.5*(omega22/omega11);
//                 const scalar bStar = (1.0/3.0)*(5.0*omega12 - omega13)/omega11;
//                 const scalar cStar = (1.0/3.0)*omega12/omega11;
                const scalar aStar = omega22/omega11;
                const scalar bStar = (5.0*omega12 - 4.0*omega13)/omega11;
                const scalar cStar = omega12/omega11;
                const scalar theta =
                    (15.0/2.0)
                   *(
                       (
                           (2.0*aStar + 5.0)
                          *(6.0*cStar - 5.0)
                          *(Wi - Wj)
                       )
                      /(
                          aStar
                         *(16.0*aStar - 12.0*bStar + 55.0)
                         *(Wi + Wj)
                       )   
                    );
                gsl_matrix_set(X, k, 0, 1.0);
                gsl_matrix_set(X, k, 1, T0);
                gsl_matrix_set(X, k, 2, T0*T0);
                gsl_matrix_set(X, k, 3, T0*T0*T0);
                gsl_vector_set(z1, k, theta);
            }
            gsl_multifit_linear(X, z1, e, cov, &chisqTherm, work);
            dictionary thermophoreticCoeffsDict;
            thermophoreticCoeffsDict.add("a", gsl_vector_get(e, 0));
            thermophoreticCoeffsDict.add("b", gsl_vector_get(e, 1));
            thermophoreticCoeffsDict.add("c", gsl_vector_get(e, 2));
            thermophoreticCoeffsDict.add("d", gsl_vector_get(e, 3));
            thermophoreticCoeffsDict.add("chiSquare", chisqTherm);
            
            thermophoreticI.add(species[j], thermophoreticCoeffsDict);
            
            scalarField DTCoeffs(4);
            forAll(DTCoeffs, s)
            {
                DTCoeffs[s] = gsl_vector_get(e, s);
            }
            Polynomial<4> DTPoly(DTCoeffs);
            word DTName = species[i] + "." + species[j];
            for (label k = 0; k < step; k++)
            {
                DTTemp[j][k] = DTPoly.value(T[k]);
            }
            
            makeGraph(T, DTTemp[j], DTName, DTPath, gFormat);
        }
        thermophoreticPropertiesDict.add(species[i], thermophoreticI);
    }    
    
    
    fileName muPath(runTime.rootPath()/runTime.caseName()/"transport"/"viscosity");
    fileName kappaPath
    (
        runTime.rootPath()/runTime.caseName()/"transport"/"conductivity" 
    );
    mkDir(muPath);
    mkDir(kappaPath);
    forAll(species, i)
    {
        makeGraph(T, mu[i], species[i], muPath, gFormat);
        makeGraph(T, kappa[i], species[i], kappaPath, gFormat);
    }        
    
    transportData.Foam::regIOobject::write();
    viscosityPropertiesDict.Foam::regIOobject::write();
    conductivityPropertiesDict.Foam::regIOobject::write();
    diffusivityPropertiesDict.Foam::regIOobject::write();
    thermophoreticPropertiesDict.Foam::regIOobject::write();
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
