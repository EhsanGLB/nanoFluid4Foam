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

#include "uniformVelocityWaveform4NFFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniformVelocityWaveform4NFFoam::uniformVelocityWaveform4NFFoam
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    amplitude_(0),
    frequency_(0),
    lenDir_(1, 0, 0),
    fourierSeriesType_("real"),
    cosCoeffs_(0),
    sinCoeffs_(0),
    real_(0),
    imaginary_(0)
{}


uniformVelocityWaveform4NFFoam::uniformVelocityWaveform4NFFoam
(
    const uniformVelocityWaveform4NFFoam& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    lenDir_(ptf.lenDir_),
    fourierSeriesType_(ptf.fourierSeriesType_),
    cosCoeffs_(ptf.cosCoeffs_),
    sinCoeffs_(ptf.sinCoeffs_),
    real_(ptf.real_),
    imaginary_(ptf.imaginary_)
{}


uniformVelocityWaveform4NFFoam::uniformVelocityWaveform4NFFoam
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    frequency_(readScalar(dict.lookup("frequency"))),
    lenDir_(dict.lookup("lenDir")),
    fourierSeriesType_(dict.lookup("fourierSeriesType")),
    cosCoeffs_(),
    sinCoeffs_(),
    real_(),
    imaginary_()
{

    if (fourierSeriesType_ == "real")
    {
        Istream& iscosCoeffs_ = dict.lookup("cosCoeffs");
        iscosCoeffs_.format(IOstream::ASCII);
        iscosCoeffs_ >> cosCoeffs_;

        Istream& issinCoeffs_ = dict.lookup("sinCoeffs");
        issinCoeffs_.format(IOstream::ASCII);
        issinCoeffs_ >> sinCoeffs_;
    }
    else if (fourierSeriesType_ == "complex")
    {
        Istream& isreal_ = dict.lookup("real");
        isreal_.format(IOstream::ASCII);
        isreal_ >> real_;

        Istream& isimaginary_ = dict.lookup("imaginary");
        isimaginary_.format(IOstream::ASCII);
        isimaginary_ >> imaginary_;
    }
    else
    {
        Info << "Select real or complex for fourierSeriesType." << endl;
        FatalErrorIn("Error") << exit(FatalError);
    }

    lenDir_ /= mag(lenDir_);

    evaluate();
}


uniformVelocityWaveform4NFFoam::uniformVelocityWaveform4NFFoam
(
    const uniformVelocityWaveform4NFFoam& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    amplitude_(fcvpvf.amplitude_),
    frequency_(fcvpvf.frequency_),
    lenDir_(fcvpvf.lenDir_),
    fourierSeriesType_(fcvpvf.fourierSeriesType_),
    cosCoeffs_(fcvpvf.cosCoeffs_),
    sinCoeffs_(fcvpvf.sinCoeffs_),
    real_(fcvpvf.real_),
    imaginary_(fcvpvf.imaginary_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::uniformVelocityWaveform4NFFoam::velocityRealFunction(scalar t)
{
    scalar omega_ = 2.0*mathematicalConstant::pi*frequency_;
    int n_ = 1;
    scalar Vt_ = cosCoeffs_[0] + sinCoeffs_[0];

    while( n_ < cosCoeffs_.size() )
    {
        Vt_ += cosCoeffs_[n_]*cos(n_*omega_*t) + sinCoeffs_[n_]*sin(n_*omega_*t);
        n_ += 1;
    }

    return amplitude_*Vt_;
}


Foam::complex Foam::uniformVelocityWaveform4NFFoam::velocityComplexFunction(scalar t)
{
    scalar omega_ = 2.0*mathematicalConstant::pi*frequency_;
    int n_ = 1;
    complex Vt_(real_[0], imaginary_[0]);
    complex i_(0.0, 1.0);

    while( n_ < real_.size() )
    {
        Vt_ += complex(real_[n_], imaginary_[n_])*exp(i_*n_*omega_*t);
        n_ += 1;
    }

    return amplitude_*Vt_;
}


Foam::scalar Foam::uniformVelocityWaveform4NFFoam::velocityFunction(scalar t)
{
    scalar f_ = 0.0;
    if (fourierSeriesType_ == "real")
    {
        f_ = velocityRealFunction(t);
    }

    else if (fourierSeriesType_ == "complex")
    {
        f_ = velocityComplexFunction(t).Re();
    }
 
    else
    {
        f_ = 0.0;
    }
    return f_;
}


void uniformVelocityWaveform4NFFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar t_ = this->db().time().timeOutputValue();

    vectorField::operator=(lenDir_*velocityFunction(t_));
    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void uniformVelocityWaveform4NFFoam::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("amplitude") << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("lenDir") << lenDir_ << token::END_STATEMENT << nl;

    os.writeKeyword("fourierSeriesType") << fourierSeriesType_ << token::END_STATEMENT << nl;
    if (fourierSeriesType_ == "real")
    {
        os.writeKeyword("cosCoeffs") << cosCoeffs_ << token::END_STATEMENT << nl;
        os.writeKeyword("sinCoeffs") << sinCoeffs_ << token::END_STATEMENT << nl;
    }

    else if (fourierSeriesType_ == "complex")
    {
        os.writeKeyword("real") << real_ << token::END_STATEMENT << nl;
        os.writeKeyword("imaginary") << imaginary_ << token::END_STATEMENT << nl;
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, uniformVelocityWaveform4NFFoam);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
