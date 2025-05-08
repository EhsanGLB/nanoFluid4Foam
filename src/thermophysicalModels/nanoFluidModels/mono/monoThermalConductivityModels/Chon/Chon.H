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
    Foam::monoThermalConductivityModels

Description
    The Chon model for mono thermal conductivity.

SourceFiles
    Chon.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef Chon_H
#define Chon_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "monoThermalConductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{

/*---------------------------------------------------------------------------*\
                           Class Chon Declaration
\*---------------------------------------------------------------------------*/

class Chon
:
    public monoThermalConductivityModel
{
    // Private data
        dictionary ChonCoeffs_;
        dimensionedScalar lbf_;
        dimensionedScalar dbf_;
        dimensionedScalar dnp_;


public:

    //- Runtime type information
    TypeName("Chon");


    // Constructors

        //- construct from components
        Chon
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T,
            const volScalarField& alpha,
            const autoPtr<baseFluid>& baseFluidPtr,
            const autoPtr<particleModel>& particleModelPtr
        );


    // Destructor

        virtual ~Chon()
        {}


    // Member Functions

        //- Return the thermal conductivity
        virtual const volScalarField  kappa() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
