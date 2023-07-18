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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::molecularTransport::molecularTransport()
:
fTrans_(1.0),
fRot_(1.0),
fVib_(1.0),
Acond_(1.0),
Bcond_(1.0),
Zrot_(1.0),
A11(1.06036),
B11(0.15610),
C11(0.19300),
D11(0.47635),
E11(1.03587),
F11(1.52996),
G11(1.76474),
H11(3.89411),
R11(0.0),
S11(0.0),
W11(0.0),
P11(0.0),
A12(1.00220),
B12(0.15530),
C12(0.16105),
D12(0.72751),
E12(0.86125),
F12(2.06848),
G12(1.95162),
H12(4.84492),
R12(0.0),
S12(0.0),
W12(0.0),
P12(0.0),
A13(0.96573),
B13(0.15611),
C13(0.44067),
D13(1.52420),
E13(2.38981),
F13(5.08063),
G13(0.0),
H13(0.0),
R13(-5.373e-4),
S13(19.2866),
W13(-1.30775),
P13(6.58711), 
A22(1.16145),
B22(0.14874),
C22(0.52487),
D22(0.77320),
E22(2.16178),
F22(2.43787),
G22(0.0),
H22(0.0),
R22(-6.435e-4),
S22(18.0323),
W22(-0.76830),
P22(7.27371)
{}

    
    
