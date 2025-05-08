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
    Foam::thermophoresisDiffusivityModels

Description
    The Brock model for thermophoresis diffusivity.

SourceFiles
    Brock.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef Brock_H
#define Brock_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "thermophoresisDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophoresisDiffusivityModels
{

/*---------------------------------------------------------------------------*\
                           Class Brock Declaration
\*---------------------------------------------------------------------------*/

class Brock
:
    public thermophoresisDiffusivityModel
{
    // Private data
        dictionary BrockCoeffs_;
        scalar Cs_;
        scalar Cm_;
        scalar Ct_;
        scalar lbf_;
        scalarField dnp_;
        PtrList<volScalarField> DTs_;


public:

    //- Runtime type information
    TypeName("Brock");


    // Constructors

        //- construct from components
        Brock
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T,
            const PtrList<volScalarField>& alphasPtr,
            const autoPtr<baseFluid>& baseFluidPtr,
            const PtrList<particleModel>& particlesProperties
        );


    // Destructor

        virtual ~Brock()
        {}


    // Member Functions

        //- Return the thermal conductivity
        inline const PtrList<volScalarField>& DTs() const
        {
            return DTs_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermophoresisDiffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
