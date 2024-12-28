/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Class
    WKBC4NFFoamFvPatchVectorField

Description
    Windkessel pressure outflow boundary condition, 3 element. 

SourceFiles
    WKBC4NFFoamFvPatchVectorField.C

Authors
    Andris Piebalgs, Imperial College London, 2017  All rights reserved
    Boram Gu, Emily Manchester, Imperial College London, 2020 All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef WKBC4NFFoam_H
#define WKBC4NFFoam_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class pulsatileVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class WKBC4NFFoam
:
    public fixedValueFvPatchScalarField
{
    // Private data

    // Name of the outlet
    scalar index_;


public:

    //- Runtime type information
    TypeName("WKBC4NFFoam");


    // Constructors

        //- Construct from patch and internal field
        WKBC4NFFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        WKBC4NFFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given WKBC4NFFoam
        //  onto a new patch
        WKBC4NFFoam
        (
            const WKBC4NFFoam&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new WKBC4NFFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        WKBC4NFFoam
        (
            const WKBC4NFFoam&,
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
                new WKBC4NFFoam(*this, iF)
            );
        }


    // Member functions

     // Return outlet index
     scalar& index()
     {
       return index_;
     }


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
