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
    parabolicRadialVelocityBC4NF

Description
    Boundary condition specifies a parabolic radial velocity inlet profile
    , given mean velocity value, radius of inlet,
    and flow direction.

    @verbatim
        inlet
        {
            type                    parabolicRadialVelocityBC4NF;
            meanValue               12.3;
            radius                  0.25;
            lenDir                  (1 0 0);
        }
    @endverbatim

SourceFiles
    parabolicRadialVelocityBC4NF.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef parabolicRadialVelocityBC4NF_H
#define parabolicRadialVelocityBC4NF_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class parabolicRadialVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class parabolicRadialVelocityBC4NF
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- mean velocity magnitude
        scalar meanValue_;

        //- Flow direction
        vector lenDir_;

        //- Radius
        scalar radius_;


public:

    //- Runtime type information
    TypeName("parabolicRadialVelocityBC4NF");


    // Constructors

        //- Construct from patch and internal field
        parabolicRadialVelocityBC4NF
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicRadialVelocityBC4NF
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given parabolicRadialVelocityBC4NF
        //  onto a new patch
        parabolicRadialVelocityBC4NF
        (
            const parabolicRadialVelocityBC4NF&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicRadialVelocityBC4NF(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicRadialVelocityBC4NF
        (
            const parabolicRadialVelocityBC4NF&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicRadialVelocityBC4NF(*this, iF)
            );
        }


    // Member functions

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
