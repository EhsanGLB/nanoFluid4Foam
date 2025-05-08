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
    The Krieger5 model for mono dynamic viscosity.

SourceFiles
    Krieger5.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef Krieger5_H
#define Krieger5_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "baseFluidDynamicViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace baseFluidDynamicViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Krieger5 Declaration
\*---------------------------------------------------------------------------*/

class Krieger5
:
    public baseFluidDynamicViscosityModel
{
    // Private data
        dictionary Krieger5Coeffs_;
        dimensionedScalar a_;
        dimensionedScalar b_;
        dimensionedScalar c_;
        dimensionedScalar beta_;
        dimensionedScalar lambda_;
        dimensionedScalar nuK_;
        dimensionedScalar muPlasma_;
        dimensionedScalar Hcrit_;


public:

    //- Runtime type information
    TypeName("Krieger5");


    // Constructors

        //- construct from components
        Krieger5
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T
        );


    // Destructor

        virtual ~Krieger5()
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
