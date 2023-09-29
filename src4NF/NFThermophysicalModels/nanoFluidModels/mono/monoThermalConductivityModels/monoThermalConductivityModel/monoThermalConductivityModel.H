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

Namespace
    Foam::monoThermalConductivityModels

Description
    A namespace for various mono thermal conductivity model implementations.

Class
    Foam::monoThermalConductivityModel

Description
    An abstract base class for incompressible monoThermalConductivityModels

SourceFiles
    monoThermalConductivityModel.C
    newmonoThermalConductivityModel.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef monoThermalConductivityModel_H
#define monoThermalConductivityModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

#include "nanoFluidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class monoThermalConductivityModel Declaration
\*---------------------------------------------------------------------------*/

class monoThermalConductivityModel
{

protected:

    // Protected data
        dictionary nanoFluidPropertiesDict_;
        const volVectorField& U_;
        const volScalarField& p_;
        const volScalarField& T_;
        const volScalarField& alpha_;
        const autoPtr<baseFluid>& baseFluidPtr_;
        const autoPtr<particleModel>& particleModelPtr_;

	dimensionedScalar kB_;

    // Private Member Functions

        //- Disallow copy construct
        monoThermalConductivityModel(const monoThermalConductivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const monoThermalConductivityModel&);


public:

    //- Runtime type information
    TypeName("monoThermalConductivityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            monoThermalConductivityModel,
            dictionary,
            (
                const dictionary& nanoFluidPropertiesDict,
                const volVectorField& U,
                const volScalarField& p,
                const volScalarField& T,
                const volScalarField& alpha,
                const autoPtr<baseFluid>& baseFluidPtr,
                const autoPtr<particleModel>& particleModelPtr
            ),
            (nanoFluidPropertiesDict, U, p, T, alpha, baseFluidPtr, particleModelPtr)
        );


    // Selectors

        //- Return a reference to the selected particle model
        static autoPtr<monoThermalConductivityModel> New
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T,
            const volScalarField& alpha,
            const autoPtr<baseFluid>& baseFluidPtr,
            const autoPtr<particleModel>& particleModelPtr
        );


    // Constructors

        //- Construct from components
        monoThermalConductivityModel
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

        virtual ~monoThermalConductivityModel()
        {}


    // Member Functions

        //- Return the thermal conductivity
        virtual const volScalarField kappa() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
