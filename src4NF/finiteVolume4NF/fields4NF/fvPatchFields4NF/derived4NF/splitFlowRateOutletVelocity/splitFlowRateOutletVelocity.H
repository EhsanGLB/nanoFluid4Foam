/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::splitFlowRateOutletVelocity

Group
    grpOutletBoundaryConditions

Description
    Velocity outlet boundary condition which corrects the extrapolated velocity
    to match the flow rate of the specified corresponding inlet patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        inletPatch   | Corresponding inlet patch name | yes  |
        volumetric   | Set volumetric or mass flow-rate | no | false
        rho          | density field name      | no          | rho
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                splitFlowRateOutletVelocity;
        inletPatch          inlet;
        value               uniform (0 0 0);
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    splitFlowRateOutletVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef splitFlowRateOutletVelocity_H
#define splitFlowRateOutletVelocity_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class splitFlowRateOutletVelocity Declaration
\*---------------------------------------------------------------------------*/

class splitFlowRateOutletVelocity
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet patch name from which the corresponding flow rate is obtained
        word inletPatchName_;

        //- Is volumetric?
        bool volumetric_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Flow split
        float flowSplit_;


    // Private member functions

        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues
        (
            const label inletPatchID,
            const RhoType& rhoOutlet,
            const RhoType& rhoInlet
        );


public:

   //- Runtime type information
   TypeName("splitFlowRateOutletVelocity");


   // Constructors

        //- Construct from patch and internal field
        splitFlowRateOutletVelocity
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        splitFlowRateOutletVelocity
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  splitFlowRateOutletVelocity
        //  onto a new patch
        splitFlowRateOutletVelocity
        (
            const splitFlowRateOutletVelocity&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        splitFlowRateOutletVelocity
        (
            const splitFlowRateOutletVelocity&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new splitFlowRateOutletVelocity(*this)
            );
        }

        //- Construct as copy setting internal field reference
        splitFlowRateOutletVelocity
        (
            const splitFlowRateOutletVelocity&,
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
                new splitFlowRateOutletVelocity(*this, iF)
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
