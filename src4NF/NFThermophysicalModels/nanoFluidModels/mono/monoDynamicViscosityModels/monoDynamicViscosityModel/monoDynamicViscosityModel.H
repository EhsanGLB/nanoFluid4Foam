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
    Foam::monoDynamicViscosityModels

Description
    A namespace for various incompressible mono dynamic viscosity model implementations.

Class
    Foam::monoDynamicViscosityModel

Description
    An abstract base class for monoDynamicViscosityModels

SourceFiles
    monoDynamicViscosityModel.C
    newmonoDynamicViscosityModel.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef monoDynamicViscosityModel_H
#define monoDynamicViscosityModel_H

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
                         Class monoDynamicViscosityModel Declaration
\*---------------------------------------------------------------------------*/

class monoDynamicViscosityModel
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
        monoDynamicViscosityModel(const monoDynamicViscosityModel&);

        //- Disallow default bitwise assignment
        void operator=(const monoDynamicViscosityModel&);


public:

    //- Runtime type information
    TypeName("monoDynamicViscosityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            monoDynamicViscosityModel,
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
        static autoPtr<monoDynamicViscosityModel> New
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
        monoDynamicViscosityModel
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

        virtual ~monoDynamicViscosityModel()
        {}


    // Member Functions

        //- Return the dynamic viscosity
        virtual const volScalarField mu() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
