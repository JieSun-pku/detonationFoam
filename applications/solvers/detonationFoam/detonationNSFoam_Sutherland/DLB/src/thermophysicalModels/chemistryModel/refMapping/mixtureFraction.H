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

Class
    Foam::mixtureFraction

Description
    Calculates and outputs the mixture fraction calculated between the
    fuel and the oxidizer.
    
    \f
        Z = (beta - beta_0)/(beta_1 - beta_0)
    \f]

    with
    \f[
        beta = 2 \frac{Z_C}{W_C} + 0.5 \frac{Z_H}{W_H} - 1 \frac{Z_O}{W_O}
    \f]

    and
    \f[
        Z_k = \sum\limits_{i=1}^N(\frac{a_{k,i} W_k}{W_i} Y_i)
    \f]

    where
    \vartable
        Z       | Mixture fraction
        beta    | Conserved scalar
        beta_0  | Conserved scalar in the oxidiser stream
        beta_1  | Conserved scalar in the fuel stream
        Z_k     | Element mass fractions (k = [C, H, O])
        W_k     | Atomic weights (k = [C, H, O])
        a_{k,i} | Number of atoms of element 'k' in specie 'i'
        N       | Number of species
        Y_i     | Species mass fractions
    \endvartable

References:

    \verbatim
        Bilger, R.; Stårner, S. & Kee, R.
        On reduced mechanisms for methane-air combustion in nonpremixed flames
        Combustion and Flame, 1990, 80, 135 - 149

        J. C. Sutherland, P. J. Smith, And J. H. Chen
        Quantification of differential diffusion in nonpremixed systems
        Combustion Theory and Modelling, Vol. 9, No. 2, May 2005, 365–383

        "Large Eddy Simulation of Fuel Spray Combustion"
        Armin Wehrfritz,
        PhD thesis, 2020, Aalto University, Finland, 2016
    \endverbatim

SourceFiles
    mixtureFraction.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureFraction_H
#define mixtureFraction_H

#include "IOdictionary.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// Mixture fraction headers
#include "atomicWeights.H"
#include "psiReactionThermo.H"

namespace Foam
{

class mixtureFraction
{

public:
    mixtureFraction() = default;

    // Construct from dict and composition 
    mixtureFraction(const dictionary& mixFracDict, const basicSpecieMixture& composition);

    scalar massFractionToMixtureFraction(const scalarField& massFraction) const;


private:

    dictionary mixFracDict_;
    wordList species_;
    List<scalar> alpha_;
    List<scalar> beta_;
    
    //- Initialize the alpha and beta values needed for 
    // mixture fraction calculation
    void initialize(const basicSpecieMixture& composition);


};

} // namespace Foam

#endif