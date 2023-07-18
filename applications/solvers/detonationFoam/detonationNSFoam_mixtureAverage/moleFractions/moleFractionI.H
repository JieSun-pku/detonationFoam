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

inline Foam::PtrList<Foam::volScalarField>&
Foam::moleFraction::X()
{
    return X_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::moleFraction::X() const
{
    return X_;
}


inline Foam::volScalarField& Foam::moleFraction::X(const label i)
{
    return X_[i];
}


inline const Foam::volScalarField& Foam::moleFraction::X
(
    const label i
) const
{
    return X_[i];
}


inline Foam::volScalarField& Foam::moleFraction::X
(
    const word& specieName
)
{
    return X_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::moleFraction::X
(
    const word& specieName
) const
{
    return X_[species_[specieName]];
}

inline void Foam::moleFraction::update()
{
    sum_ *= 0.0;
    forAll(X_, i)
    {
        X_[i] = Y_[i]/W_[i];
	sum_ += X_[i];
    }
    forAll(X_, i)
    {
        X_[i] /= sum_;
    }
}
// ************************************************************************* //
