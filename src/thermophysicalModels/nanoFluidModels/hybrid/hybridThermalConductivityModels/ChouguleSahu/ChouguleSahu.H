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
    Foam::hybridThermalConductivityModel

Description
    The ChouguleSahu for hybrid thermal conductivity.

SourceFiles
    ChouguleSahu.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef ChouguleSahu_H
#define ChouguleSahu_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "hybridThermalConductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace hybridThermalConductivityModels
{

/*---------------------------------------------------------------------------*\
                           Class ChouguleSahu Declaration
\*---------------------------------------------------------------------------*/

class ChouguleSahu
:
    public hybridThermalConductivityModel
{
    // Private data
        dictionary ChouguleSahuCoeffs_;
        dimensionedScalar dbf_;
        scalarList dnpList_;


public:

    //- Runtime type information
    TypeName("ChouguleSahu");


    // Constructors

        //- construct from components
        ChouguleSahu
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

        virtual ~ChouguleSahu()
        {}


    // Member Functions

        //- Return the thermal conductivity
        virtual const volScalarField kappa() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hybridThermalConductivityModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
