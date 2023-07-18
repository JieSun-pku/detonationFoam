/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/blttkgl/DLBFoam

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

#include "mixtureFraction.H"



// Constructor
Foam::mixtureFraction::mixtureFraction(const dictionary& mixFracDict, const basicSpecieMixture& composition)
    : mixFracDict_(mixFracDict), species_(composition.species()), alpha_(composition.species().size(), 0.0),
      beta_(2, 0.0) 
{
    initialize(composition);
}

void Foam::mixtureFraction::initialize(const basicSpecieMixture& composition)
{
    forAll(alpha_, i)
    {
        const dictionary& dict =
            mixFracDict_.subDict(species_[i]).subDict("elements");
        scalar a0(
            2.0 * dict.lookupOrDefault<label>("C", 0) / composition.W(i));
        scalar a1(
            0.5 * dict.lookupOrDefault<label>("H", 0) / composition.W(i));
        scalar a2(
            -1.0 * dict.lookupOrDefault<label>("O", 0) / composition.W(i));
        alpha_[i] = a0 + a1 + a2;
    }

    List<List<scalar>> Yconst(2, List<scalar>(species_.size(), 0.0));
    forAll(species_, i)
    {
        Yconst[0][i] = mixFracDict_.subDict("oxidizerMassFractions")
                           .lookupOrDefault<scalar>(species_[i], 0.0);
        Yconst[1][i] = mixFracDict_.subDict("fuelMassFractions")
                           .lookupOrDefault<scalar>(species_[i], 0.0);
    }

    scalar YoxTot = 0.0;
    scalar YfuTot = 0.0;
    forAll(species_, i)
    {
        YoxTot += Yconst[0][i];
        YfuTot += Yconst[1][i];
    }

    if(mag(1.0 - YoxTot) > SMALL || mag(1.0 - YfuTot) > SMALL)
    {
        FatalErrorIn("createMixtureFraction.H :")
            << "oxidizerMassFractions or fuelMassFractions do not sum up to 1.0"
            << abort(FatalError);
    }

    forAll(species_, i)
    {
        beta_[0] += alpha_[i] * Yconst[0][i]; // oxidizer
        beta_[1] += alpha_[i] * Yconst[1][i]; // fuel
    }


    // Stoichiometric mixture fraction
    scalar Z_st = (0.0 - beta_[0])/(beta_[1] - beta_[0]);

    Info<<"Stoichiomentric mixture fraction Zst = " << Z_st << endl;

}



Foam::scalar Foam::mixtureFraction::massFractionToMixtureFraction(const scalarField& massFraction) const
{

    scalar beta = 0.0; 
    forAll(massFraction, iField)
    {
        beta += alpha_[iField] * massFraction[iField];
    }

    return (beta - beta_[0]) / (beta_[1] - beta_[0]);



}
