/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::ROUNDLimiter

Description
    Class with limiter function which returns the limiter for the
    ROUND differencing scheme based on phict obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    ROUND.C

\*---------------------------------------------------------------------------*/

#ifndef ROUNDL_H
#define ROUNDL_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ROUNDLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class ROUNDLLimiter
:
    public LimiterFunc
{

public:

    ROUNDLLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
  

       scalar r = LimiterFunc::r
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );      
 
          scalar rr=0.0;
          
         if(r>0.0){
          
          scalar a=(r)/(1.0+r);
          scalar a8=sqr(sqr(sqr(1.0+5.0*sqr(-1.0+a))));
          scalar a4=sqr(sqr(1.0+12.0*sqr(a)));
          
          scalar maxal=max(18000.0*sqr(sqr(0.97-a))*(0.97-a)*sqr(-0.55+a)*(-0.55+a),0.0);
          scalar maxbl=max(1100.0*sqr(0.47-a)*(0.47-a)*sqr(-0.05+a)*(-0.05+a),0.0);
          scalar maxsum=1.0/3.0+5.0*a/6.0+maxal+maxbl;
      
          rr=2.0*(-a+(0.5+0.5*a)/a8+(1.0-1.0/a8)*(3.0*a/2.0/a4+(1.0-1.0/a4)*maxsum))*(1.0+r);

        }
      
      
       return rr;
        
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
