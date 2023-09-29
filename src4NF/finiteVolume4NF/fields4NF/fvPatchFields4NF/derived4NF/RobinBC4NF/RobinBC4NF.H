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
    Foam::RobinBC4NF

Description
    Boundary condition specifies a parabolic radial velocity inlet profile
    , given mean velocity value, radius of inlet, flow direction.

    @verbatim
        inlet
        {
            type                    RobinBC4NF;
            flux                    100.0;
            ho                      0.1;
        }
    @endverbatim

SourceFiles
    RobinBC4NF.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef RobinBC4NF_H
#define RobinBC4NF_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformRobinFvPatch Declaration
\*---------------------------------------------------------------------------*/

class RobinBC4NF
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Name of the thermal conductivity field
        word kappaEffName_;

        //- Flux
        scalarList flux_;

        //- direction
        vector dir_;

        //- convective heat transfer
        scalar ho_;

        //- Out temperature
        scalar To_;


public:

    //- Runtime type information
    TypeName("RobinBC4NF");


    // Constructors

        //- Construct from patch and internal field
        RobinBC4NF
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        RobinBC4NF
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given RobinBC4NF onto
        //  a new patch
        RobinBC4NF
        (
            const RobinBC4NF&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        RobinBC4NF
        (
            const RobinBC4NF&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new RobinBC4NF(*this)
            );
        }

        //- Construct as copy setting internal field reference
        RobinBC4NF
        (
            const RobinBC4NF&,
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
                new RobinBC4NF(*this, iF)
            );
        }


    // Member functions

        //- Calculate the flux(location)
        virtual scalar fluxFunction(scalar location);

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
