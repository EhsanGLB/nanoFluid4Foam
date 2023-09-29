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
    lumpedLIBBC4NF

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    lumpedLIBBC4NF.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef lumpedLIBBC4NF_H
#define lumpedLIBBC4NF_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformTemperatureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class lumpedLIBBC4NF
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- state
        word state_;

        //- LIB density
        scalar rho_;

        //- LIB heat specific capacity
        scalar Cp_;

        //- current
        scalarList current_;

        //- resistance
        scalarList resistance_;

        //- entropy coefficient
        scalar dEdT_;

        //- Name of the temperature field
        word TName_;

        //- Name of the thermal conductivity field
        word kappaEffName_;

        //- LIB density
        scalar initTemp_;

        //- state of charge
        scalar SOC_;

        //- LIB capacity
        scalar capacity_;

        //- LIB volume
        scalar volume_;

        //- currentForm
        word currentForm_;

        //- Amplitude
        scalar amplitude_;

        //- Frequency
        scalar frequency_;

        //- Real component
        scalarList real_;

        //- Imaginary component
        scalarList imaginary_;



public:

    //- Runtime type information
    TypeName("lumpedLIBBC4NF");


    // Constructors

        //- Construct from patch and internal field
        lumpedLIBBC4NF
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        lumpedLIBBC4NF
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given lumpedLIBBC4NF
        //  onto a new patch
        lumpedLIBBC4NF
        (
            const lumpedLIBBC4NF&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new lumpedLIBBC4NF(*this)
            );
        }

        //- Construct as copy setting internal field reference
        lumpedLIBBC4NF
        (
            const lumpedLIBBC4NF&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new lumpedLIBBC4NF(*this, iF)
            );
        }

        //- Construct as copy
        lumpedLIBBC4NF
        (
            const lumpedLIBBC4NF&
        );


    // Member functions

        //- Calculate the resistance(LIBTemperature)
        virtual scalar resistanceFunction(scalar LIBTemp);

        //- Calculate the currentPolynomialFunction(time)
        virtual scalar currentPolynomialFunction(scalar t);

        //- Calculate the currentWavyFunction(time)
        virtual complex currentWavyFunction(scalar t);

        //- Calculate the currentFunction(time)
        virtual scalar currentFunction(scalar t);

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
