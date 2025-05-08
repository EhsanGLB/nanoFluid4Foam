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
    Foam::monoDynamicViscosityModels

Description
    The muPolyExp model for mono dynamic viscosity.

SourceFiles
    muPolyExp.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef muPolyExp_H
#define muPolyExp_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "monoDynamicViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace monoDynamicViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class muPolyExp Declaration
\*---------------------------------------------------------------------------*/

class muPolyExp
:
    public monoDynamicViscosityModel
{
    // Private data
        dictionary muPolyExpCoeffs_;
        scalar a0_;
        scalar a1_;
        scalar a2_;
        scalar a3_;
        scalar a4_;
        scalar a5_;
        scalar a6_;
        scalar b_;


public:

    //- Runtime type information
    TypeName("muPolyExp");


    // Constructors

        //- construct from components
        muPolyExp
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

        virtual ~muPolyExp()
        {}


    // Member Functions

        //- Return the dynamic viscosity
        virtual const volScalarField mu() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoDynamicViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
