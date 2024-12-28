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
    Foam::brownianDiffusivityModels

Description
    A namespace for various brownian diffusivity model implementations.

Class
    Foam::brownianDiffusivityModel

Description
    An abstract base class for brownianDiffusivityModels

SourceFiles
    brownianDiffusivityModel.C
    newbrownianDiffusivityModel.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef brownianDiffusivityModel_H
#define brownianDiffusivityModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

#include "baseFluid.H"//-nanoFluid4Foam
#include "particleModel.H"//-nanoFluid4Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class brownianDiffusivityModel Declaration
\*---------------------------------------------------------------------------*/

class brownianDiffusivityModel
{

protected:

    // Protected data
        dictionary nanoFluidPropertiesDict_;
        const volVectorField& U_;
        const volScalarField& p_;
        const volScalarField& T_;
        const PtrList<volScalarField>& alphasPtr_;
        const autoPtr<baseFluid>& baseFluidPtr_;
        const PtrList<particleModel>& particlesProperties_;

	dimensionedScalar kB_;
        dimensionedScalar dd_;



    // Private Member Functions

        //- Disallow copy construct
        brownianDiffusivityModel(const brownianDiffusivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const brownianDiffusivityModel&);


public:

    //- Runtime type information
    TypeName("brownianDiffusivityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            brownianDiffusivityModel,
            dictionary,
            (
                const dictionary& nanoFluidPropertiesDict,
                const volVectorField& U,
                const volScalarField& p,
                const volScalarField& T,
                const PtrList<volScalarField>& alphasPtr,
                const autoPtr<baseFluid>& baseFluidPtr,
                const PtrList<particleModel>& particlesProperties
            ),
            (nanoFluidPropertiesDict, U, p, T, alphasPtr, baseFluidPtr, particlesProperties)
        );


    // Selectors

        //- Return a reference to the selected particle model
        static autoPtr<brownianDiffusivityModel> New
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T,
            const PtrList<volScalarField>& alphasPtr,
            const autoPtr<baseFluid>& baseFluidPtr,
            const PtrList<particleModel>& particlesProperties
        );


    // Constructors

        //- Construct from components
        brownianDiffusivityModel
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

        virtual ~brownianDiffusivityModel()
        {}


    // Member Functions

        //- Return the thermal conductivity
        virtual const PtrList<volScalarField>& DBs() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
