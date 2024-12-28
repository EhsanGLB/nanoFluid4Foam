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
    Foam::brownianDiffusivityModels

Description
    An incompressible Brownian baseFluid model.

SourceFiles
    Brownian.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef Brownian_H
#define Brownian_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "brownianDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace brownianDiffusivityModels
{

/*---------------------------------------------------------------------------*\
                           Class Brownian Declaration
\*---------------------------------------------------------------------------*/

class Brownian
:
    public brownianDiffusivityModel
{
    // Private data
        dictionary BrownianCoeffs_;
        scalarField dnp_;
        volScalarField mubf_;
        PtrList<volScalarField> DBs_;


public:

    //- Runtime type information
    TypeName("Brownian");


    // Constructors

        //- construct from components
        Brownian
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

        virtual ~Brownian()
        {}


    // Member Functions

        //- Return the thermal conductivity
        inline const PtrList<volScalarField>& DBs() const
        {
            return DBs_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace brownianDiffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
