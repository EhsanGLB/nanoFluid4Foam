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

Class
    Foam::newMaterialConst

Description
    newMaterialConst

SourceFiles
    newMaterialConst.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef newMaterialConst_H
#define newMaterialConst_H

#include "material.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace materials
{

/*---------------------------------------------------------------------------*\
                           Class newMaterialConst Declaration
\*---------------------------------------------------------------------------*/

class newMaterialConst
:
    public material
{
    // Private data

        scalar rho_;
        scalar K_;
        scalar cp_;
        scalar beta_;
        scalar rhoR_;
        scalar mu_;
        scalar pv_;


public:

    //- Runtime type information
    TypeName("newMaterialConst");


    // Constructors

        //- Construct from Istream
        newMaterialConst(Istream& is);


    // Member Functions

            //- material density [kg/m^3]
            inline scalar rho(scalar T) const;

            //- material thermal conductivity  [W/(m K)]
            inline scalar K(scalar T) const;

            //- material heat capacity [J/(kg K)]
            inline scalar cp(scalar T) const;

            //- material thermal expansion [1/K]
            inline scalar beta(scalar T) const;

            //- material electrical resistivity [ohm.m]
            inline scalar rhoR(scalar T) const;

            //- material dynamic viscosity [Pa s]
            inline scalar mu(scalar T) const;

            //- Vapour pressure [Pa]
            inline scalar pv(scalar T) const;
};

/*
    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            material::writeData(os); os << nl;
            rho.writeData(os); os << nl;
            os  << rho_ << token::SPACE
                << K_ << token::SPACE
                << cp_ << token::SPACE
                << beta_ << token::SPACE
                << rhoR_ << token::SPACE
                << mu_ << token::SPACE
                << pv_;
        }


        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const newMaterialConst& m)
        {
            m.writeData(os);
            return os;
        }
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "newMaterialConstI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
