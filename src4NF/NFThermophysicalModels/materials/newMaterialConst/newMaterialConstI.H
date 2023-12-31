/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

namespace Foam
{
namespace materials
{

inline scalar newMaterialConst::rho(scalar T) const
{
    return rho_;
}


inline scalar newMaterialConst::K(scalar T) const
{
    return K_;
}


inline scalar newMaterialConst::cp(scalar T) const
{
    return cp_;
}


inline scalar newMaterialConst::beta(scalar T) const
{
    return beta_;
}


inline scalar newMaterialConst::rhoR(scalar T) const
{
    return rhoR_;
}


inline scalar newMaterialConst::mu(scalar T) const
{
    return mu_;
}


inline scalar newMaterialConst::pv(scalar T) const
{
    return pv_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materials
} // End namespace Foam

// ************************************************************************* //
