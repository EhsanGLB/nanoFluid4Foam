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
    Foam::convCoeffModels

Description
    An incompressible flatPlateFC nanoFluid model.

SourceFiles
    flatPlateFC.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef flatPlateFC_H
#define flatPlateFC_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "convCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace convCoeffModels
{

/*---------------------------------------------------------------------------*\
                           Class flatPlateFC Declaration
\*---------------------------------------------------------------------------*/

class flatPlateFC
:
    public convCoeffModel
{
        //- Out velocity
        scalar Vo_;

        //- Length
        scalar L_;



public:

    //- Runtime type information
    TypeName("flatPlateFC");


    // Constructors

        //- construct from components
        flatPlateFC
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~flatPlateFC()
        {}


    // Member Functions

        //- Reynolds number [-]
        virtual scalar Re() const;

        //- Return the average convective heat transfer coefficient [W/m^2 K]
        virtual scalar hoAvg(scalar avgT) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace convCoeffModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
