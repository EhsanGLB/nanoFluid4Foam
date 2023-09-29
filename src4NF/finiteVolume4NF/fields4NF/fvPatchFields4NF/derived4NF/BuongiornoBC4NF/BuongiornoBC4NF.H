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
    Foam::BuongiornoBC4NF

Description
    Foam::BuongiornoBC4NF

SourceFiles
    BuongiornoBC4NF.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef BuongiornoBC4NFs_H
#define BuongiornoBC4NFs_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformBuongiornoFvPatch Declaration
\*---------------------------------------------------------------------------*/

class BuongiornoBC4NF
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Name of the alpha field
        word alphaName_;

        //- Name of the alpha field
        word particleName_;

        //- massFlux
        scalar massFlux_;

        //- convective heat transfer
        scalar ho_;

        //- Out temperature
        scalar alphao_;


public:

    //- Runtime type information
    TypeName("BuongiornoBC4NF");


    // Construcalphaors

        //- Construct from patch and internal field
        BuongiornoBC4NF
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        BuongiornoBC4NF
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BuongiornoBC4NF onalphao
        //  a new patch
        BuongiornoBC4NF
        (
            const BuongiornoBC4NF&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        BuongiornoBC4NF
        (
            const BuongiornoBC4NF&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new BuongiornoBC4NF(*this)
            );
        }

        //- Construct as copy setting internal field reference
        BuongiornoBC4NF
        (
            const BuongiornoBC4NF&,
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
                new BuongiornoBC4NF(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
