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
    Foam::mixtureFractionRefMapper

Description
    Extends the base class RefMapperBase by implementing a
    reference mapping method using Bilger's mixture fraction definition,
    calculated between fuel and oxidizer.

References
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
    mixtureFractionRefMapper.C

See also
    mixtureFraction.H
    
\*---------------------------------------------------------------------------*/

#ifndef mixtureFractionRefMapper_H
#define mixtureFractionRefMapper_H

#include "IOdictionary.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "psiReactionThermo.H"

#include "mixtureFraction.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class mixtureFractionRefMapper 
{

public:
    mixtureFractionRefMapper() = default;

    mixtureFractionRefMapper(
        const dictionary& dict, const basicSpecieMixture& composition)
        : dict_(dict), coeffsDict_(dict.subDict("refmapping")),
          active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
          Ztolerance_(coeffsDict_.lookupOrDefault<scalar>("tolerance", 1e-4)),
          Ttolerance_(coeffsDict_.lookupOrDefault<scalar>("deltaT", VGREAT)),
          mixture_fraction_()
    {
        if (active())
        {
            mixture_fraction_ = mixtureFraction(
              coeffsDict_.subDict("mixtureFractionProperties"),
              composition);
        }
    }



    //- Check if the massFraction is within Ztolerance_
    bool shouldMap(const scalarField& massFraction) const;


    //- Check if the temperature is within Ttolerance_
    bool temperatureWithinRange(scalar Ti, scalar Tref) const;


    //- Is reference mapping active?
    bool active() const
    {
        return active_;
    }


private:
    const dictionary dict_;
    const dictionary coeffsDict_;

    // Is reference mapping active?
    Switch active_;

    // Tolerance for reference mapping
    scalar Ztolerance_;

    // Tolerance for temperature checking (in Kelvins)
    scalar Ttolerance_;
    
    // Mixture fraction object
    mixtureFraction mixture_fraction_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
