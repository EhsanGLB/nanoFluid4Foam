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
    Foam::baseFluidDynamicViscosityModels

Description
    The Yeleswarapu model for mono dynamic viscosity.

SourceFiles
    Yeleswarapu.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef Yeleswarapu_H
#define Yeleswarapu_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "baseFluidDynamicViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace baseFluidDynamicViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Yeleswarapu Declaration
\*---------------------------------------------------------------------------*/

class Yeleswarapu
:
    public baseFluidDynamicViscosityModel
{
    // Private data
        dictionary YeleswarapuCoeffs_;
        dimensionedScalar a1_;
        dimensionedScalar a2_;
        dimensionedScalar a3_;
        dimensionedScalar b1_;
        dimensionedScalar b2_;
        dimensionedScalar b3_;
        dimensionedScalar k_;
        dimensionedScalar muPlasma_;


public:

    //- Runtime type information
    TypeName("Yeleswarapu");


    // Constructors

        //- construct from components
        Yeleswarapu
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T
        );


    // Destructor

        virtual ~Yeleswarapu()
        {}


    // Member Functions

        //- Return the dynamic viscosity
        virtual const volScalarField mu() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace baseFluidDynamicViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
