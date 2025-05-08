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
    Foam::robin4NFFoam

Description
    Boundary condition specifies a parabolic radial velocity inlet profile
    , given mean velocity value, radius of inlet, flow direction.

    @verbatim
        hotWall
        {
            type                    robin4NFFoam;
            Kappa                   kappa;
            qo                      (q0 q1 q2);
            ho                      (h0 h1 h2);
            To                      (T0 T1 T2);
        }
    @endverbatim

SourceFiles
    robin4NFFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef robin4NFFoam_H
#define robin4NFFoam_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformRobinFvPatch Declaration
\*---------------------------------------------------------------------------*/

class robin4NFFoam
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Name of the thermal conductivity field
        word kappaName_;

        //- Flux depends on time
        scalarList qo_;

        //- convective heat transfer
        scalarList ho_;

        //- Out temperature
        scalarList To_;


public:

    //- Runtime type information
    TypeName("robin4NFFoam");


    // Constructors

        //- Construct from patch and internal field
        robin4NFFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        robin4NFFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given robin4NFFoam onto
        //  a new patch
        robin4NFFoam
        (
            const robin4NFFoam&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        robin4NFFoam
        (
            const robin4NFFoam&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new robin4NFFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        robin4NFFoam
        (
            const robin4NFFoam&,
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
                new robin4NFFoam(*this, iF)
            );
        }


    // Member functions

        //- Calculate the qo(time)
        virtual scalar qoFunction(scalar time);

        //- Calculate the ho(time)
        virtual scalar hoFunction(scalar time);

        //- Calculate the To(time)
        virtual scalar ToFunction(scalar time);

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
