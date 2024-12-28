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
    Foam::material

Description
    The thermophysical properties of a material

SourceFiles
    material.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef material_H
#define material_H

#include "scalar.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

#include "colors.H"//-nanoFluid4Foam
#include "expPolyfunc.H"//-nanoFluid4Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class material Declaration
\*---------------------------------------------------------------------------*/

class material
{
    // Private data

        //- Molecular weight [kg/kmol]
        scalar W_;

        //- Critical temperature [K]
        scalar Tc_;

        //- Critical pressure [Pa]
        scalar Pc_;

        //- Critical volume [m^3/mol]
        scalar Vc_;

        //- Melting point [K]
        scalar Tm_;

        //- Melting range [K]
        scalar Tmr_;

        //- Melting enthalpy [J/kg]
        scalar hm_;


public:

    TypeName("material");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            material,
            ,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            material,
            Istream,
            (Istream& is),
            (is)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            material,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

    // Constructors

        //- Construct from components
        material(){}

        //material(Istream& is){}
        material
        (
            scalar W,
            scalar Tc,
            scalar Pc,
            scalar Vc,
            scalar Tm,
            scalar Tmr,
            scalar hm
        )
        :
            W_(W),
            Tc_(Tc),
            Pc_(Pc),
            Vc_(Vc),
            Tm_(Tm),
            Tmr_(Tmr),
            hm_(hm)
        {}

        //- Construct from Istream
        material(Istream& is)
        :
            W_(readScalar(is)),
            Tc_(readScalar(is)),
            Pc_(readScalar(is)),
            Vc_(readScalar(is)),
            Tm_(readScalar(is)),
            Tmr_(readScalar(is)),
            hm_(readScalar(is))
        {}

        //- Return a pointer to a new material created from Istream
        static autoPtr<material> New(Istream& is);

        //- Return a pointer to a new material created from dictionary
        static autoPtr<material> New(const dictionary& dict);


    //- Destructor
    virtual ~material()
    {}


    // Member Functions

        // Physical constants
            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- Critical temperature [K]
            inline scalar Tc() const;

            //- Critical pressure [Pa]
            inline scalar Pc() const;

            //- Critical volume [m^3/mol]
            inline scalar Vc() const;

            //- Melting point [K]
            inline scalar Tm() const;

            //- Melting range [K]
            inline scalar Tmr() const;

            //- Melting enthalpy [J/kg]
            inline scalar hm() const;

            //- Melting entropy [J/kg.K]
            inline scalar sm() const;

        // Physical property pure virtual functions

            //- material density [kg/m^3]
            virtual scalar rho(scalar T) const = 0;

            //- material thermal conductivity [W/(m K)]
            virtual scalar K(scalar T) const = 0;

            //- material heat capacity [J/(kg K)]
            virtual scalar cp(scalar T) const = 0;

            //- material thermal expansion [1/K]
            virtual scalar beta(scalar T) const = 0;

            //- material electrical resistivity [ohm.m]
            virtual scalar rhoR(scalar T) const = 0;

            //- material dynamic viscosity [Pa s]
            virtual scalar mu(scalar T) const = 0;

            //- Vapour pressure [Pa]
            virtual scalar pv(scalar T) const = 0;


    // I-O

        //- Write the function coefficients
        virtual void writeData(Ostream& os) const
        {
            os  << W_ << token::SPACE
                << Tc_ << token::SPACE
                << Pc_ << token::SPACE
                << Vc_ << token::SPACE
                << Tm_ << token::SPACE
                << Tmr_ << token::SPACE
                << hm_;
        }

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const material& m)
        {
            m.writeData(os);
            return os;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "materialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
