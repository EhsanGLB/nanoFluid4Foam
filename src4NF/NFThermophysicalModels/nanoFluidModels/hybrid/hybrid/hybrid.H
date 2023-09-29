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
    An incompressible hybrid nanoFluid model.

SourceFiles
    hybrid.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef hybrid_H
#define hybrid_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "nanoFluidModel.H"
#include "hybridThermalConductivityModel.H"
#include "hybridDynamicViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nanoFluidModels
{

/*---------------------------------------------------------------------------*\
                           Class hybrid Declaration
\*---------------------------------------------------------------------------*/

class hybrid
:
    public nanoFluidModel
{

        //volScalarField alpha_;
        autoPtr<hybridThermalConductivityModel> hybridThermalConductivityModelPtr_;
        autoPtr<hybridDynamicViscosityModel> hybridDynamicViscosityModelPtr_;


public:

    //- Runtime type information
    TypeName("hybrid");

    // Constructors

        //- construct from components
        hybrid
        (
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T,
            const PtrList<volScalarField>& alphasPtr
        );


    // Destructor

        virtual ~hybrid()
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

        //- Return the total alpha
        virtual const volScalarField alpha() const;

        //- print
        virtual void print() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nanoFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
