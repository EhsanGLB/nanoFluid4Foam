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
    Foam::particleModels

Description
    A namespace for various incompressible particleModel implementations.

Class
    Foam::particleModel

Description
    An abstract base class for incompressible particleModels

SourceFiles
    particleModel.C
    newparticleModel.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef particleModel_H
#define particleModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "volFields.H"

#include "materialMixture.H"//-nanoFluid4Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class particleModel Declaration
\*---------------------------------------------------------------------------*/

class particleModel
{

protected:

    // Protected data
        wordList particleList_;
        dictionary nanoFluidPropertiesDict_;
        const volVectorField& U_;
        const volScalarField& p_;
        const volScalarField& T_;

        word particleName_;
        dictionary particleDict_;
        dictionary particleSubDict_;


    // Private Member Functions

        //- Disallow copy construct
        particleModel(const particleModel&);

        //- Disallow default bitwise assignment
        void operator=(const particleModel&);


public:

    //- Runtime type information
    TypeName("particleModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            particleModel,
            dictionary,
            (
                const wordList& particleList,
                const dictionary& nanoFluidPropertiesDict,
                const volVectorField& U,
                const volScalarField& p,
                const volScalarField& T

            ),
            (particleList, nanoFluidPropertiesDict, U, p, T)
        );


    // Selectors

        //- Return a reference to the selected particle model
        static autoPtr<particleModel> New
        (
            const wordList& particleList,
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T
        );


    // Constructors


        //- Construct from components
        particleModel
        (
            const wordList& particleList,
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T
        );


    // Destructor

        virtual ~particleModel()
        {}


    // Member Functions

        //- Return the particleList
        const wordList particleList() const
        {
            return particleList_;
        }

        //- Return the nanoFluidPropertiesDict
        const dictionary nanoFluidPropertiesDict() const
        {
            return nanoFluidPropertiesDict_;
        }

        //- Return the U
        const volVectorField U() const
        {
            return U_;
        }

        //- Return the p
        const volScalarField p() const
        {
            return p_;
        }

        //- Return the T
        const volScalarField T() const
        {
            return T_;
        }

        //- Return the name
        const word particleName() const
        {
            return particleName_;
        }

    // Member Functions

        //- Return the density [kg/m^3]
        virtual const volScalarField rho() const = 0;

        //- Return the thermal conductivity [W/(m K)]
        virtual const volScalarField kappa() const = 0;

        //- Return the specific heat capacity [J/(kg K)]
        virtual const volScalarField Cp() const = 0;

        //- Return the thermal expansion [1/K]
        virtual const volScalarField beta() const = 0;

        //- Return the electrical resistivity [ohm.m]
        virtual const volScalarField rhoR() const = 0;

        //- Return the dynamic viscosity [Pa s]
        virtual const volScalarField mu() const = 0;

        //- Return the vapour pressure [Pa]
        virtual const volScalarField pv() const = 0;

        //- Return the melting enthalpy [J/kg]
        virtual const volScalarField hm() const = 0;

        //- Return the melting point [K]
        virtual const volScalarField Tm() const = 0;

        //- Return the melting range [K]
        virtual const volScalarField Tmr() const = 0;

        //- Correct the molten
        virtual const volScalarField molten() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
