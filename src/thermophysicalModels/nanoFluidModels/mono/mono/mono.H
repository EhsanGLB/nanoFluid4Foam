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
    Foam::nanoFluidModels

Description
    An incompressible mono nanoFluid model.

SourceFiles
    mono.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef mono_H
#define mono_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "nanoFluidModel.H"
#include "monoThermalConductivityModel.H"
#include "monoDynamicViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nanoFluidModels
{

/*---------------------------------------------------------------------------*\
                           Class mono Declaration
\*---------------------------------------------------------------------------*/

class mono
:
    public nanoFluidModel
{

        autoPtr<particleModel> particleModelPtr_;
        autoPtr<monoThermalConductivityModel> monoThermalConductivityModelPtr_;
        autoPtr<monoDynamicViscosityModel> monoDynamicViscosityModelPtr_;


public:

    //- Runtime type information
    TypeName("mono");


    // Constructors

        //- construct from components
        mono
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T,
            const PtrList<volScalarField>& alphasPtr
        );


    // Destructor

        virtual ~mono()
        {}


    // Member Functions

        //- Return the density [kg/m^3]
        virtual const volScalarField rho() const;

        //- Return the thermal conductivity [W/(m K)]
        virtual const volScalarField kappa() const;

        //- Return the specific heat capacity [J/(kg K)]
        virtual const volScalarField Cp() const;

        //- Return the thermal expansion [1/K]
        virtual const volScalarField beta() const;

        //- Return the electrical resistivity [ohm.m]
        virtual const volScalarField rhoR() const;

        //- Return the dynamic viscosity [Pa s]
        virtual const volScalarField mu() const;

        //- Return the vapour pressure [Pa]
        virtual const volScalarField pv() const;

        //- Return the melting enthalpy [J/kg]
        virtual const volScalarField hm() const;

        //- Return the melting point [K]
        virtual const volScalarField Tm() const;

        //- Return the melting range [K]
        virtual const volScalarField Tmr() const;

        //- print
        virtual void print() const;

        //- Return the alpha
        virtual const volScalarField alpha() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nanoFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
