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
    Foam::nanoFluidModels

Description
    A namespace for various incompressible nanoFluidModel implementations.

Class
    Foam::nanoFluidModel

Description
    An abstract base class for incompressible nanoFluidModels

SourceFiles
    nanoFluidModel.C
    newnanoFluidModel.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef nanoFluidModel_H
#define nanoFluidModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

#include "baseFluid.H"//-nanoFluid4Foam
#include "particleModel.H"//-nanoFluid4Foam
#include "brownianDiffusivityModel.H"//-nanoFluid4Foam
#include "thermophoresisDiffusivityModel.H"//-nanoFluid4Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class nanoFluidModel Declaration
\*---------------------------------------------------------------------------*/

class nanoFluidModel
{
protected:

    // Protected data
        dictionary nanoFluidPropertiesDict_;
        const volVectorField& U_;
        const volScalarField& p_;
        const volScalarField& T_;
        const PtrList<volScalarField>& alphasPtr_;

        wordList nanoFluidModelList_;
        word nanoFluidModelType_;

        autoPtr<baseFluid> baseFluidPtr_;

        //- The names of the particles
        List<wordList> particlesComponents_;

        //- The particle properties
        PtrList<particleModel> particlesProperties_;


        PtrList<volScalarField> DBs_;
        //autoPtr<brownianDiffusivityModel> brownianDiffusivityModelPtr_;

        PtrList<volScalarField> DTs_;
        //autoPtr<thermophoresisDiffusivityModel> thermophoresisDiffusivityModelPtr_;



    // Private Member Functions

        //- Disallow copy construct
        nanoFluidModel(const nanoFluidModel&);

        //- Disallow default bitwise assignment
        void operator=(const nanoFluidModel&);


public:

    //- Runtime type information
    TypeName("nanoFluidModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            nanoFluidModel,
            dictionary,
            (
                const dictionary& nanoFluidPropertiesDict,
                const volVectorField& U,
                const volScalarField& p,
                const volScalarField& T,
                const PtrList<volScalarField>& alphasPtr
            ),
            (nanoFluidPropertiesDict, U, p, T, alphasPtr)
        );


    // Selectors

        //- Return a reference to the selected baseFluid model
        static autoPtr<nanoFluidModel> New
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T,
            const PtrList<volScalarField>& alphasPtr
        );


    // Constructors

        //- Construct from components
        nanoFluidModel
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T,
            const PtrList<volScalarField>& alphasPtr
        );


    // Destructor

        virtual ~nanoFluidModel()
        {}


    // Member Functions

        //- Return the nanoFluidModelType
        bool nanoFluidModelType()
        {
            return ((nanoFluidModelType_ == "Homogenous") ? false : true);
        }

        //- Return the baseFluid
        inline baseFluid& baseFluidPtr()
        {
            return baseFluidPtr_();
        }

        //- Return the particles
        inline const PtrList<particleModel>& particlesProperties() const
        {
            return particlesProperties_;
        }

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

        //- Read nanoFluidPropertiesDict dictionary
        virtual bool read(const dictionary& nanoFluidPropertiesDict);

        //- print
        virtual void print() const = 0;

        //- Return the Brownian Diffusivity [m^2/s]
        inline const PtrList<volScalarField>& DBs() const
        {
            return DBs_;
        }

        //- Return the thermophoresis Diffusivity [m^2/s]
        virtual const PtrList<volScalarField>& DTs() const
        {
            return DTs_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
