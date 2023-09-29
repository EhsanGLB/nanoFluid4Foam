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
    Foam::newMaterialFuncT

Description
    newMaterialFuncT

SourceFiles
    newMaterialFuncT.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef newMaterialFuncT_H
#define newMaterialFuncT_H

#include "material.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materials
{

/*---------------------------------------------------------------------------*\
                           Class newMaterialFuncT Declaration
\*---------------------------------------------------------------------------*/

class newMaterialFuncT
:
    public material
{
    // Private data
        dictionary dict_;
        NFExpPolyFunc rho_;
        NFExpPolyFunc K_;
        NFExpPolyFunc cp_;
        NFExpPolyFunc beta_;
        NFExpPolyFunc rhoR_;
        NFExpPolyFunc mu_;
        NFExpPolyFunc pv_;


public:

    //- Runtime type information
    TypeName("newMaterialFuncT");


    // Constructors

        //- Construct from Istream
        newMaterialFuncT(const dictionary& dict);


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
            rho_.writeData(os); os << nl;
            K_.writeData(os); os << nl;
            cp_.writeData(os); os << nl;
            beta_.writeData(os); os << nl;
            rhoR_.writeData(os); os << nl;
            mu_.writeData(os); os << nl;
            pv_.writeData(os); os << nl;
        }


        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const newMaterialFuncT& m)
        {
            m.writeData(os);
            return os;
        }
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "newMaterialFuncTI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
