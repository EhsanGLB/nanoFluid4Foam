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
    Foam::materialMixture

Description
    A mixture of materials.

    The theory in this class is to a very large extent based on the book:
    'The Properties of Gases & materials' 4th ed.
    by Reid, Prausnitz and Poling

    For now it does not do much, since the perfect gas equation is used.

    The dictionary constructor searches for the entry @c materialComponents,
    which is a wordList. The material properties of each component can either
    be contained within a @c materialProperties sub-dictionary or (for legacy
    purposes) can be found directly in the dictionary.
    The @c materialProperties sub-dictionary entry should be used when possible
    to avoid conflicts with identically named gas-phase entries.

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef materialMixture_H
#define materialMixture_H

#include "word.H"
#include "scalarField.H"
#include "PtrList.H"
#include "autoPtr.H"

#include "material.H"//-nanoFluid4Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class materialMixture Declaration
\*---------------------------------------------------------------------------*/

class materialMixture
{
    // Private data

        // maximum reduced temperature
        static const scalar TrMax;

        //- The names of the materials
        wordList components_;

        //- The material properties
        PtrList<material> properties_;


public:


    // Constructors

        //- Construct from dictionary
        materialMixture
        (
            const wordList& components,
            const dictionary& propertiesSubDict
        );


    // Selectors

        //- Select construct from dictionary
        static autoPtr<materialMixture> New
        (
            const wordList& components,
            const dictionary& propertiesSubDict
        );


    // Member Functions

        inline const wordList& components() const
        {
            return components_;
        }

        //- Return the material properties
        inline const PtrList<material>& properties() const
        {
            return properties_;
        }

        //- Calculate the critical temperature of mixture
        scalar Tc(const scalarField& x) const;

        //- Calculate the mean molecular weight [kg/kmol]
        //  from mole fractions
        scalar W(const scalarField& x) const;

        //- Returns the mass fractions, given mole fractions
        scalarField Y(const scalarField& X) const;

        //- Returns the mole fractions, given mass fractions
        scalarField X(const scalarField& Y) const;

        //- Calculate the mixture density [kg/m^3]
        scalar rho
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Estimate thermal conductivity  [W/(m K)]
        // Li's method, Eq. 10-12.27 - 10.12-19
        scalar K
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture heat capacity [J/(kg K)]
        scalar cp
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture thermal expansion [1/K]
        scalar beta
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture electrical resistivity [ohm.m]
        scalar rhoR
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture dynamic viscosity [Pa s]
        scalar mu
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture vapour pressure [Pa]
        scalar pv
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture melting enthalpy [J/kg]
        scalar hm
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture melting entropy [J/kg.K]
        scalar sm
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture melting point [K]
        scalar Tm
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture melting range [K]
        scalar Tmr
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
