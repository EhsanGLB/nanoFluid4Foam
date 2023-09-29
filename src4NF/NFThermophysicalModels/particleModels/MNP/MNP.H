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
    Foam::particleModels::MNP

Description
    An incompressible MNP(Micro/Nano Particle) model.

SourceFiles
    MNP.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef MNP_H
#define MNP_H

#include "dimensionedScalar.H"
#include "volFields.H"

#include "particleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace particleModels
{

/*---------------------------------------------------------------------------*\
                           Class MNP Declaration
\*---------------------------------------------------------------------------*/

class MNP
:
    public particleModel
{
    // Private data
        wordList MNPComponents_;
        scalarField MNPMoleFraction_;
        autoPtr<materialMixture> materialMixturePtr_;


public:

    //- Runtime type information
    TypeName("MNP");


    // Constructors

        //- construct from components
        MNP
        (
            const wordList& particleList,
            const dictionary& nanoFluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T
        );


    // Destructor

        virtual ~MNP()
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

        //- Return the molten [-]
        virtual const volScalarField molten() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace particleModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
