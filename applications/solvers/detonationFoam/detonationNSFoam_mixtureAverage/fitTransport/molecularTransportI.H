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
 #include "molecularTransport.H"


//---------------
inline Foam::scalar Foam::molecularTransport::F(scalar LJ,const scalar T)
{
    
    scalar eByT_= LJ/T;
    scalar eBy298_= LJ/298.0;
    const scalar pi_ = Foam::constant::mathematical::pi;
    scalar F =
    (
        1.0 + 0.5*Foam::pow(pi_,1.5)*Foam::sqrt(eBy298_)
       +eBy298_*(0.25*pi_*pi_+2.0)
       +Foam::pow(pi_*eBy298_,1.5)
    )
   /(
	1.0+0.5*Foam::pow(pi_,1.5)*Foam::sqrt(eByT_)
       +eByT_*(0.25*pi_*pi_+2.0)
       +Foam::pow(pi_*eByT_,1.5)
    );
    return F;	
 }
//-Collision Integrals
inline Foam::scalar Foam::molecularTransport::CollisionIntegral11
(
    const scalar Treduced
)
{
    
    scalar omega11;
    omega11 =
        A11/Foam::pow(Treduced,B11)
       +C11*Foam::exp(-min(600.0, D11*Treduced))
       +E11*Foam::exp(-min(600.0, F11*Treduced))
       +G11*Foam::exp(-min(600.0, H11*Treduced));
    return omega11;
 
}
inline Foam::scalar Foam::molecularTransport::CollisionIntegral22
(
    const scalar Treduced
)
{
   scalar omega22;
   omega22 =
       A22/Foam::pow(Treduced,B22)
      +C22*Foam::exp(-min(600.0, D22*Treduced))
      +E22*Foam::exp(-min(600.0, F22*Treduced))
      +R22*Foam::pow(Treduced,B22)*Foam::sin(S22*Foam::pow(Treduced,W22)-P22);
   return omega22;
}

inline Foam::scalar Foam::molecularTransport::CollisionIntegral12
(
    const scalar Treduced
)
{
   scalar omega12;
   omega12 =
       A12/Foam::pow(Treduced,B12)
      +C12*Foam::exp(-min(600.0, D12*Treduced))
      +E12*Foam::exp(-min(600.0, F12*Treduced))
      +G12*Foam::exp(-min(600.0, H12*Treduced));
   return omega12;
}

inline Foam::scalar Foam::molecularTransport::CollisionIntegral13
(
    const scalar Treduced
)
{
   scalar omega13;
   omega13 =
       A13/Foam::pow(Treduced,B13)
      +C13*Foam::exp(-min(600.0, D13*Treduced))
      +E13*Foam::exp(-min(600.0, F13*Treduced))
      +R13*Foam::pow(Treduced,B13)*Foam::sin(S13*Foam::pow(Treduced,W13)-P13);
   return omega13;
}
inline Foam::scalar Foam::molecularTransport::viscosity
(
    const scalar m,
    const scalar T,
    const scalar sigma,
    const scalar omega22
)
{
    scalar visc;
    visc = 
        5.0/16.0
       *(
           Foam::sqrt
           (
               Foam::constant::mathematical::pi
              *m
              *Foam::constant::physicoChemical::k.value()
              *T
           )
	  /(Foam::constant::mathematical::pi*sigma*sigma*omega22)
        );
  return visc;    
}
inline Foam::scalar Foam::molecularTransport::diffusivity
(
    const scalar m,
    const scalar T,
    const scalar sigma,
    const scalar omega11,
    const scalar p
)
{
    scalar diff;
    diff =
        3.0/16.0
       *(
           Foam::sqrt
           (
               2.0*Foam::constant::mathematical::pi
              *Foam::pow((Foam::constant::physicoChemical::k.value()*T),3.0)
              /m               
           )
	  /(Foam::constant::mathematical::pi*sigma*sigma*omega11*p)
	);
    return diff;	     
}

//*************** //
