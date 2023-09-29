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

#include "parabolicTroughCollectorBC4NF.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "scalar.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicTroughCollectorBC4NF::parabolicTroughCollectorBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    kappaEffName_("kappaEff"),
    range_(0),
    LCR1_(0),
    LCR2_(0),
    LCR3_(0),
    glassRadProp_(0),
    tubeRadProp_(0),
    glassInDiam_(0.0),
    glassOutDiam_(0.0),
    tubeInDiam_(0.0),
    tubeOutDiam_(0.0),
    PTCLength_(0.0),
    longDir_(0, 0, 0),
    radDir_(0, 0, 0),
    ambConv_(0.0),
    ambTemp_(0.0),
    IntRad_(0.0)
{}


Foam::parabolicTroughCollectorBC4NF::parabolicTroughCollectorBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaEffName_(dict.lookupOrDefault<word>("kappaEff", "kappaEff")),
    range_(),
    LCR1_(),
    LCR2_(),
    LCR3_(),
    glassRadProp_(),
    tubeRadProp_(),
    glassInDiam_(readScalar(dict.lookup("glassInDiam"))),
    glassOutDiam_(readScalar(dict.lookup("glassOutDiam"))),
    tubeInDiam_(readScalar(dict.lookup("tubeInDiam"))),
    tubeOutDiam_(readScalar(dict.lookup("tubeOutDiam"))),
    PTCLength_(readScalar(dict.lookup("PTCLength"))),
    longDir_(dict.lookup("longDir")),
    radDir_(dict.lookup("radDir")),
    ambConv_(readScalar(dict.lookup("ambConv"))),
    ambTemp_(readScalar(dict.lookup("ambTemp"))),
    IntRad_(readScalar(dict.lookup("IntRad")))
{
        Istream& isrange_ = dict.lookup("range");
        isrange_.format(IOstream::ASCII);
        isrange_ >> range_;

        Istream& isLCR1_ = dict.lookup("LCR1");
        isLCR1_.format(IOstream::ASCII);
        isLCR1_ >> LCR1_;

        Istream& isLCR2_ = dict.lookup("LCR2");
        isLCR2_.format(IOstream::ASCII);
        isLCR2_ >> LCR2_;

        Istream& isLCR3_ = dict.lookup("LCR3");
        isLCR3_.format(IOstream::ASCII);
        isLCR3_ >> LCR3_;

        Istream& isglassRadProp_ = dict.lookup("glassRadProp");
        isglassRadProp_.format(IOstream::ASCII);
        isglassRadProp_ >> glassRadProp_;

        Istream& istubeRadProp_= dict.lookup("tubeRadProp");
        istubeRadProp_.format(IOstream::ASCII);
        istubeRadProp_ >> tubeRadProp_;

    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::parabolicTroughCollectorBC4NF::parabolicTroughCollectorBC4NF
(
    const parabolicTroughCollectorBC4NF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    range_(ptf.range_),
    LCR1_(ptf.LCR1_),
    LCR2_(ptf.LCR2_),
    LCR3_(ptf.LCR3_),
    glassRadProp_(ptf.glassRadProp_),
    tubeRadProp_(ptf.tubeRadProp_),
    glassInDiam_(ptf.glassInDiam_),
    glassOutDiam_(ptf.glassOutDiam_),
    tubeInDiam_(ptf.tubeInDiam_),
    tubeOutDiam_(ptf.tubeOutDiam_),
    PTCLength_(ptf.PTCLength_),
    longDir_(ptf.longDir_),
    radDir_(ptf.radDir_),
    ambConv_(ptf.ambConv_),
    ambTemp_(ptf.ambTemp_),
    IntRad_(ptf.IntRad_)
{}


Foam::parabolicTroughCollectorBC4NF::parabolicTroughCollectorBC4NF
(
    const parabolicTroughCollectorBC4NF& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    range_(ptf.range_),
    LCR1_(ptf.LCR1_),
    LCR2_(ptf.LCR2_),
    LCR3_(ptf.LCR3_),
    glassRadProp_(ptf.glassRadProp_),
    tubeRadProp_(ptf.tubeRadProp_),
    glassInDiam_(ptf.glassInDiam_),
    glassOutDiam_(ptf.glassOutDiam_),
    tubeInDiam_(ptf.tubeInDiam_),
    tubeOutDiam_(ptf.tubeOutDiam_),
    PTCLength_(ptf.PTCLength_),
    longDir_(ptf.longDir_),
    radDir_(ptf.radDir_),
    ambConv_(ptf.ambConv_),
    ambTemp_(ptf.ambTemp_),
    IntRad_(ptf.IntRad_)
{}


Foam::parabolicTroughCollectorBC4NF::parabolicTroughCollectorBC4NF
(
    const parabolicTroughCollectorBC4NF& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    range_(ptf.range_),
    LCR1_(ptf.LCR1_),
    LCR2_(ptf.LCR2_),
    LCR3_(ptf.LCR3_),
    glassRadProp_(ptf.glassRadProp_),
    tubeRadProp_(ptf.tubeRadProp_),
    glassInDiam_(ptf.glassInDiam_),
    glassOutDiam_(ptf.glassOutDiam_),
    tubeInDiam_(ptf.tubeInDiam_),
    tubeOutDiam_(ptf.tubeOutDiam_),
    PTCLength_(ptf.PTCLength_),
    longDir_(ptf.longDir_),
    radDir_(ptf.radDir_),
    ambConv_(ptf.ambConv_),
    ambTemp_(ptf.ambTemp_),
    IntRad_(ptf.IntRad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::parabolicTroughCollectorBC4NF::LCRFunction(scalar teta)
{
////range coefficients
    scalar teta0_(range_[0]);
    scalar teta1_(range_[1]);
    scalar teta2_(range_[2]);
    scalar teta3_(range_[3]);

    scalar f_ = 0.0;

    if ( teta0_ <= teta && teta < teta1_)
    {
        int n1_ = 0;

        while( n1_ < LCR1_.size() )
        {
            f_ += LCR1_[n1_]*pow(teta, n1_);
            n1_ += 1;
        }
    }

    else if ( teta1_ <= teta && teta < teta2_)
    {
        int n2_ = 0;

        while( n2_ < LCR2_.size() )
        {
            f_ += LCR2_[n2_]*pow(teta, n2_);
            n2_ += 1;
        }
    }

    else if ( teta2_ <= teta && teta < teta3_)
    {
        int n3_ = 0;

        while( n3_ < LCR3_.size() )
        {
            f_ += LCR3_[n3_]*pow(teta, n3_);
            n3_ += 1;
        }
    }

    return f_;
}


void Foam::parabolicTroughCollectorBC4NF::updateCoeffs()
{
    if (updated())
    {
        return;
    }

////radiative properties of glass
    scalar glassAbs_(glassRadProp_[0]);
    scalar glassEmis_(glassRadProp_[1]);
    scalar glassTrans_(glassRadProp_[2]);

////radiative properties of tube
    scalar tubeAbs_(tubeRadProp_[0]);
    scalar tubeEmis_(tubeRadProp_[1]);
    //scalar tubeTrans_(tubeRadProp_[2]);

////constant parameters
    const scalar sigmaSB_(5.67e-8);

////fixed values
    //scalar glassInArea_( mathematicalConstant::pi * PTCLength_ * glassInDiam_ );
    scalar glassOutArea_( mathematicalConstant::pi * PTCLength_ * glassOutDiam_ );
    scalar tubeInArea_( mathematicalConstant::pi * PTCLength_ * tubeInDiam_ );
    scalar tubeOutArea_( mathematicalConstant::pi * PTCLength_ * tubeOutDiam_ );
    scalar skyTemp_( 0.0522 * pow(ambTemp_, 1.5) );

//// To calculate degree field
    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    boundBox bb_(patch().patch().localPoints(), true);
    vector ctr_ = 0.5*(bb_.max() + bb_.min());
    const vectorField& patchCent_ = patch().Cf();
    scalarField degree_(patchCent_.size(), 0.0);
    scalar t_ = -1.1*PTCLength_/2.0;
    scalar tMax_ = 1.1*PTCLength_/2.0;
    scalar tStep_ = PTCLength_/5000.0;
    scalar error_ = 0.1;

    while( t_ <= tMax_ )
    {
        vector rC_ = ctr_ + t_ * (longDir_/mag(longDir_));

        forAll(patchCent_, i)
        {
            vector rPrime_ = patchCent_[i] - rC_;
            scalar residual_ = mag( (rPrime_/mag(rPrime_)) & (longDir_/mag(longDir_)) );
            if ( (residual_ <= error_) )
            {
                degree_[i] = acos( (rPrime_/mag(rPrime_)) & (radDir_/mag(radDir_)) ) * (180/mathematicalConstant::pi);
            }
        }

        t_ += tStep_;
    }


//// To calculate LCR average and tube average temperature
    scalar LCRAvg_(0.0);
    for(int i=0; i<180; i++)
    {
        LCRAvg_ += LCRFunction(scalar(i));
    }
    LCRAvg_ /= 180;

    scalar tubeAvgTemp_(0.0);

    forAll(Tp, i)
    {
        tubeAvgTemp_ += Tp[i] * patch().magSf()[i];
    }

    tubeAvgTemp_ /= tubeInArea_;

//// To calculate heat loss from tube
    scalar x1 = 10;
    scalar x2 = 1000;
    scalar x3;
    int count = 0;
    int iter = 100;

    for (count=0; count<iter; count++)
    {
        x3 = (x1 + x2)/2;
        scalar f1 = sigmaSB_*tubeOutArea_*(pow(tubeAvgTemp_,4) - pow(x1,4)) / ( 1/tubeEmis_ + ((1-glassEmis_)/glassEmis_)*(tubeOutDiam_/glassInDiam_) ) + (IntRad_*LCRAvg_*glassAbs_*glassOutArea_) - (ambConv_*glassOutArea_*(x1-ambTemp_) + glassEmis_*sigmaSB_*glassOutArea_*(pow(x1,4) - pow(skyTemp_,4)));

        scalar f3 = sigmaSB_*tubeOutArea_*(pow(tubeAvgTemp_,4) - pow(x3,4)) / ( 1/tubeEmis_ + ((1-glassEmis_)/glassEmis_)*(tubeOutDiam_/glassInDiam_) ) + (IntRad_*LCRAvg_*glassAbs_*glassOutArea_) - (ambConv_*glassOutArea_*(x3-ambTemp_) + glassEmis_*sigmaSB_*glassOutArea_*(pow(x3,4) - pow(skyTemp_,4)));
        if(neg(f1*f3)){x2 = x3;}
        else if((mag(x1 - x2) < 1e-2)||(f3 == 1e-5)){break;}
        else{x1 = x3;}
        count++;
    }

    scalar glassAvgTemp_ = x3;
    scalar tubeGlassHeatTrans_ = sigmaSB_*tubeOutArea_*(pow(tubeAvgTemp_,4) - pow(glassAvgTemp_,4)) / ( 1/tubeEmis_ + ((1-glassEmis_)/glassEmis_)*(tubeOutDiam_/glassInDiam_) );

//// To calculate gradient
    const fvPatchScalarField& kappaEffp = lookupPatchField<volScalarField, scalar>(kappaEffName_);
    forAll(degree_, i)
    {
        gradient()[i] = ( glassTrans_*tubeAbs_*IntRad_*LCRFunction(degree_[i]) - tubeGlassHeatTrans_/tubeOutArea_)/kappaEffp[i];
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::parabolicTroughCollectorBC4NF::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("range") << range_ << token::END_STATEMENT << nl;
    os.writeKeyword("LCR1") << LCR1_ << token::END_STATEMENT << nl;
    os.writeKeyword("LCR2") << LCR2_ << token::END_STATEMENT << nl;
    os.writeKeyword("LCR3") << LCR3_ << token::END_STATEMENT << nl;
    os.writeKeyword("glassRadProp") << glassRadProp_ << token::END_STATEMENT << nl;
    os.writeKeyword("tubeRadProp") << tubeRadProp_ << token::END_STATEMENT << nl;
    os.writeKeyword("glassInDiam") << glassInDiam_ << token::END_STATEMENT << nl;
    os.writeKeyword("glassOutDiam") << glassOutDiam_ << token::END_STATEMENT << nl;
    os.writeKeyword("tubeInDiam") << tubeInDiam_ << token::END_STATEMENT << nl;
    os.writeKeyword("tubeOutDiam") << tubeOutDiam_ << token::END_STATEMENT << nl;
    os.writeKeyword("PTCLength") << PTCLength_ << token::END_STATEMENT << nl;
    os.writeKeyword("longDir") << longDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("radDir") << radDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("ambConv") << ambConv_ << token::END_STATEMENT << nl;
    os.writeKeyword("ambTemp") << ambTemp_ << token::END_STATEMENT << nl;
    os.writeKeyword("IntRad") << IntRad_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        parabolicTroughCollectorBC4NF
    );
}

// ************************************************************************* //
