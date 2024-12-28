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
    parabolicRadialVelocityWaveform4NFFoam

Description
    Boundary condition specifies a power law velocity inlet profile
    (fixed value), given mean velocity value,
    flow direction n and direction of the power law coordinate y.

SourceFiles
    parabolicRadialVelocityWaveform4NFFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef parabolicRadialVelocityWaveform4NFFoam_H
#define parabolicRadialVelocityWaveform4NFFoam_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "complex.H"
#include "complex4NFFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class parabolicRadialVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class parabolicRadialVelocityWaveform4NFFoam
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Amplitude
        scalar amplitude_;

        //- Frequency
        scalar frequency_;

        //- Flow direction
        vector longDir_;

        //- R
        scalar R_;

        //- fourierSeriesType
        word fourierSeriesType_;

        //- cosCoeffs component
        scalarList cosCoeffs_;

        //- sinCoeffs component
        scalarList sinCoeffs_;

        //- Real component
        scalarList real_;

        //- Imaginary component
        scalarList imaginary_;



public:

    //- Runtime type information
    TypeName("parabolicRadialVelocityWaveform4NFFoam");


    // Constructors

        //- Construct from patch and internal field
        parabolicRadialVelocityWaveform4NFFoam
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicRadialVelocityWaveform4NFFoam
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given parabolicRadialVelocityWaveform4NFFoam
        //  onto a new patch
        parabolicRadialVelocityWaveform4NFFoam
        (
            const parabolicRadialVelocityWaveform4NFFoam&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicRadialVelocityWaveform4NFFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicRadialVelocityWaveform4NFFoam
        (
            const parabolicRadialVelocityWaveform4NFFoam&,
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
                new parabolicRadialVelocityWaveform4NFFoam(*this, iF)
            );
        }


    // Member functions

        //- Calculate the velocityReal(time)
        virtual scalar velocityRealFunction(scalar t);

        //- Calculate the velocityComplex(time)
        virtual complex velocityComplexFunction(scalar t);

        //- Calculate the velocity(time)
        virtual scalar velocityFunction(scalar t);

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
