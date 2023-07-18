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
    Foam::molecularTransport

Description
    Calculates the transport propertiers for fitting.

SourceFiles
    molecularTransport.C
    molecularTransportI.H

\*---------------------------------------------------------------------------*/

#ifndef molecularTransport_H
#define molecularTransport_H
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class molecularTransport Declaration
\*---------------------------------------------------------------------------*/

class molecularTransport
{

protected:

    // Protected data

	scalar fTrans_;
	scalar fRot_;
	scalar fVib_;
	scalar Acond_;
	scalar Bcond_;
	scalar Zrot_;	
	scalar A11,B11,C11,D11,E11,F11,G11,H11,R11,S11,P11,W11;//CI11
	scalar A12,B12,C12,D12,E12,F12,G12,H12,R12,S12,P12,W12;//CI12
	scalar A13,B13,C13,D13,E13,F13,G13,H13,R13,S13,P13,W13;//CI13
        scalar A22,B22,C22,D22,E22,F22,G22,H22,R22,S22,P22,W22;//CI22
        scalar ratioA_;
	scalar ratioB_;
	scalar ratioC_;

        

public:

    // Constructor
        molecularTransport();
    //- Destructor
    virtual ~molecularTransport()
    {}


    // Member functions
	//-Calculate collision integrals
	inline scalar CollisionIntegral22(const scalar Treduced);
	inline scalar CollisionIntegral11(const scalar Treduced);
	inline scalar CollisionIntegral12(const scalar Treduced);
	inline scalar CollisionIntegral13(const scalar Treduced);
	inline scalar F(scalar LJ,const scalar T);
	inline scalar viscosity
	(
            const scalar m,
            const scalar T,
            const scalar sigma,
            const scalar omega22
        );
	inline scalar diffusivity
	(
            const scalar m,
            const scalar T,
            const scalar sigma,
            const scalar omega11,
            const scalar p=1e5
        );
	

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "molecularTransportI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
