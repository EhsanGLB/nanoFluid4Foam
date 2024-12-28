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

#include "parabolicRadialVelocityWaveform4NFFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parabolicRadialVelocityWaveform4NFFoam::parabolicRadialVelocityWaveform4NFFoam
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    amplitude_(0),
    frequency_(0),
    longDir_(1, 0, 0),
    R_(0),
    fourierSeriesType_("real"),
    cosCoeffs_(0),
    sinCoeffs_(0),
    real_(0),
    imaginary_(0)
{}


parabolicRadialVelocityWaveform4NFFoam::parabolicRadialVelocityWaveform4NFFoam
(
    const parabolicRadialVelocityWaveform4NFFoam& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    longDir_(ptf.longDir_),
    R_(ptf.R_),
    fourierSeriesType_(ptf.fourierSeriesType_),
    cosCoeffs_(ptf.cosCoeffs_),
    sinCoeffs_(ptf.sinCoeffs_),
    real_(ptf.real_),
    imaginary_(ptf.imaginary_)
{}


parabolicRadialVelocityWaveform4NFFoam::parabolicRadialVelocityWaveform4NFFoam
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    frequency_(readScalar(dict.lookup("frequency"))),
    longDir_(dict.lookup("longDir")),
    R_(readScalar(dict.lookup("R"))),
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

    longDir_ /= mag(longDir_);

    evaluate();
}


parabolicRadialVelocityWaveform4NFFoam::parabolicRadialVelocityWaveform4NFFoam
(
    const parabolicRadialVelocityWaveform4NFFoam& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    amplitude_(fcvpvf.amplitude_),
    frequency_(fcvpvf.frequency_),
    longDir_(fcvpvf.longDir_),
    R_(fcvpvf.R_),
    fourierSeriesType_(fcvpvf.fourierSeriesType_),
    cosCoeffs_(fcvpvf.cosCoeffs_),
    sinCoeffs_(fcvpvf.sinCoeffs_),
    real_(fcvpvf.real_),
    imaginary_(fcvpvf.imaginary_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::parabolicRadialVelocityWaveform4NFFoam::velocityRealFunction(scalar t)
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


Foam::complex Foam::parabolicRadialVelocityWaveform4NFFoam::velocityComplexFunction(scalar t)
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


Foam::scalar Foam::parabolicRadialVelocityWaveform4NFFoam::velocityFunction(scalar t)
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


void parabolicRadialVelocityWaveform4NFFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar t_ = this->db().time().timeOutputValue();
    boundBox bb_(patch().patch().localPoints(), true);
    vector ctr_ = 0.5*(bb_.max() + bb_.min());
    const vectorField& c_ = patch().Cf();
    scalarField rp_ = mag(c_ - ctr_);

    vectorField::operator=(longDir_*2.0*velocityFunction(t_)*(1.0 - pow(rp_/R_, 2.0)));
    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void parabolicRadialVelocityWaveform4NFFoam::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("amplitude") << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("longDir") << longDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;

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

makePatchTypeField(fvPatchVectorField, parabolicRadialVelocityWaveform4NFFoam);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
