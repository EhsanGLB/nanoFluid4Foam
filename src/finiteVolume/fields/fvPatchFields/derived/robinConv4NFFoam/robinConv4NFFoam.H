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
    Foam::robinConv4NFFoam

Description
    Boundary condition specifies a parabolic radial velocity inlet profile
    , given mean velocity value, radius of inlet, flow direction.

    @verbatim
        hotWall
        {
            type                    robinConv4NFFoam;
            fluxLoc                 (a0 a1 a2);
            fluxTime                (b0 b1 b2);
            ho                      0.1;
            To                      300;
        }
    @endverbatim

SourceFiles
    robinConv4NFFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef robinConv4NFFoam_H
#define robinConv4NFFoam_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

#include "convCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformRobinFvPatch Declaration
\*---------------------------------------------------------------------------*/

class robinConv4NFFoam
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        dictionary dictc_;

        //- Name of the temperature field
        word TName_;

        //- Name of the thermal conductivity field
        word kappaEffName_;

        //- Flux depends on location
        scalarList fluxLoc_;

        //- Flux depends on time
        scalarList fluxTime_;

        //- direction
        vector dir_;

        //- Out temperature
        scalar To_;

        //- convective model
        autoPtr<convCoeffModel> convCoeffModelPtr_;

public:

    //- Runtime type information
    TypeName("robinConv4NFFoam");


    // Constructors

        //- Construct from patch and internal field
        robinConv4NFFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        robinConv4NFFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given robinConv4NFFoam onto
        //  a new patch
        robinConv4NFFoam
        (
            const robinConv4NFFoam&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        robinConv4NFFoam
        (
            const robinConv4NFFoam&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new robinConv4NFFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        robinConv4NFFoam
        (
            const robinConv4NFFoam&,
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
                new robinConv4NFFoam(*this, iF)
            );
        }


    // Member functions

        //- Return the dictBC
        const dictionary& dictBC() const
        {
            return dictc_;
        }

        //- Calculate the fluxLoc(location)
        virtual scalar fluxLocFunction(scalar location);

        //- Calculate the heat transfer convection
        virtual scalar hoAvg();

        //- Calculate the fluxTime(time)
        virtual scalar fluxTimeFunction(scalar time);

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
