/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::ROUNDAplusLimiter

Description
    Class with limiter function which returns the limiter for the
    ROUNDAplus differencing scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    ROUNDAplus.C

\*---------------------------------------------------------------------------*/

#ifndef ROUNDAplus_H
#define ROUNDAplus_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ROUNDAplusLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class ROUNDAplusLimiter
:
    public LimiterFunc
{

public:

    ROUNDAplusLimiter(Istream&)
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


        scalar rr = 0.0;
        scalar beta1 = 1100;
        scalar beta2 = 800;
        scalar lambd = 0.15;
        
        if(r>1.0){

        scalar A = sqr(sqr(sqr(1.0+r)));
        scalar C = beta2*sqr(sqr(1.0-r))+32.0*sqr(sqr(1.0+r));
        scalar D = beta2*sqr(sqr(1.0-r))+16.0*sqr(sqr(1.0+r));
        
        rr=min((6.0*beta2*(1.0-lambd)*sqr(sqr(-1.0+r))*C+256.0*A*(2.0+r))/(3.0*D*D),2.0-2.0*lambd);
        
        
        }
        else if (r>0 && r<=1){ 
        
        scalar A = sqr(sqr(sqr(1.0+r)));
        scalar B = beta1*sqr(sqr(0.5-0.5*r))+sqr(sqr(1.0+r));      
        
        rr=min((A*(2.0-5.0*r)+6.0*r*B*B)/(3.0*B*B),2*r);
        
        }        

        return rr;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
