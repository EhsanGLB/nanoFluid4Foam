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

\*---------------------------------------------------------------------------*/

#include "powerLawRadialVelocityBC4NF.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

powerLawRadialVelocityBC4NF::powerLawRadialVelocityBC4NF
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    meanValue_(0),
    lenDir_(1, 0, 0),
    radius_(0.0),
    m_(5.0)
{}


powerLawRadialVelocityBC4NF::powerLawRadialVelocityBC4NF
(
    const powerLawRadialVelocityBC4NF& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    meanValue_(ptf.meanValue_),
    lenDir_(ptf.lenDir_),
    radius_(ptf.radius_),
    m_(ptf.m_)
{}


powerLawRadialVelocityBC4NF::powerLawRadialVelocityBC4NF
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    meanValue_(readScalar(dict.lookup("meanValue"))),
    lenDir_(dict.lookup("lenDir")),
    radius_(readScalar(dict.lookup("radius"))),
    m_(readScalar(dict.lookup("m")))
{
    lenDir_ /= mag(lenDir_);

    evaluate();
}


powerLawRadialVelocityBC4NF::powerLawRadialVelocityBC4NF
(
    const powerLawRadialVelocityBC4NF& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    meanValue_(fcvpvf.meanValue_),
    lenDir_(fcvpvf.lenDir_),
    radius_(fcvpvf.radius_),
    m_(fcvpvf.m_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void powerLawRadialVelocityBC4NF::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get range and orientation
    boundBox bb(patch().patch().localPoints(), true);

    vector ctr = 0.5*(bb.max() + bb.min());

    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate for the parabolic profile
    scalarField rp_ = mag(c - ctr);

    scalar Cm_ = (1 + 1/m_)*(2 + 1/m_) / 2;

    vectorField::operator=(lenDir_*Cm_*meanValue_*pow( (1.0 - (rp_/radius_)), 1/m_));
    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void powerLawRadialVelocityBC4NF::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("meanValue") << meanValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("lenDir") << lenDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius") << radius_ << token::END_STATEMENT << nl;
    os.writeKeyword("m") << m_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, powerLawRadialVelocityBC4NF);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
