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

#include "windowBC4NF.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::windowBC4NF::windowBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    kappaEffName_("kappaEff"),
    qo_(0.0),
    qoRange_(0.0),
    To_(0.0),
    ho_(0.0),
    emissivity_(0.0)
{}


Foam::windowBC4NF::windowBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaEffName_(dict.lookupOrDefault<word>("kappaEff", "kappaEff")),
    qo_(),
    qoRange_(),
    To_(),
    ho_(readScalar(dict.lookup("ho"))),
    emissivity_(readScalar(dict.lookup("emissivity")))
{
    Istream& isqo_ = dict.lookup("qo");
    isqo_.format(IOstream::ASCII);
    isqo_ >> qo_;

    Istream& isqoRange_ = dict.lookup("qoRange");
    isqoRange_.format(IOstream::ASCII);
    isqoRange_ >> qoRange_;

    Istream& isTo_ = dict.lookup("To");
    isTo_.format(IOstream::ASCII);
    isTo_ >> To_;

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


Foam::windowBC4NF::windowBC4NF
(
    const windowBC4NF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    qo_(ptf.qo_),
    qoRange_(ptf.qoRange_),
    To_(ptf.To_),
    ho_(ptf.ho_),
    emissivity_(ptf.emissivity_)
{}


Foam::windowBC4NF::windowBC4NF
(
    const windowBC4NF& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    qo_(ptf.qo_),
    qoRange_(ptf.qoRange_),
    To_(ptf.To_),
    ho_(ptf.ho_),
    emissivity_(ptf.emissivity_)
{}


Foam::windowBC4NF::windowBC4NF
(
    const windowBC4NF& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    qo_(ptf.qo_),
    qoRange_(ptf.qoRange_),
    To_(ptf.To_),
    ho_(ptf.ho_),
    emissivity_(ptf.emissivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::windowBC4NF::qoFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    if( qoRange_[0] <= (t/3600) && (t/3600) <= qoRange_[1] )
    {
        while( n_ < qo_.size() )
        {
            f_ += qo_[n_]*pow(t/3600, n_);
            n_ += 1;
        }
    }
    else{f_=0.0;}

    return f_;
}


Foam::scalar Foam::windowBC4NF::ToFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < To_.size() )
    {
        f_ += To_[n_]*pow(t/3600, n_);
        n_ += 1;
    }

    return f_;
}


void Foam::windowBC4NF::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- get the currernt time
    scalar t_ = this->db().time().value();


    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& kappaEffp = lookupPatchField<volScalarField, scalar>(kappaEffName_);
    const scalar sigmaSB_(5.67e-8);
    scalar skyTemp_( 0.0522 * pow(ToFunction(t_), 1.5) );

    forAll(patch().Cf(), i)
    {
        gradient()[i] = ( qoFunction(t_) - ho_ * (Tp[i] - ToFunction(t_)) - emissivity_*sigmaSB_*(pow(ToFunction(t_), 4.0) - pow(skyTemp_, 4.0)) )/kappaEffp[i];
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
Info << "qo = " << qoFunction(t_) << endl;
Info << "To = " << ToFunction(t_) << endl;
}


void Foam::windowBC4NF::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("qo") << qo_ << token::END_STATEMENT << nl;
    os.writeKeyword("qoRange") << qoRange_ << token::END_STATEMENT << nl;
    os.writeKeyword("To") << To_ << token::END_STATEMENT << nl;
    os.writeKeyword("ho") << ho_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivity") << emissivity_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        windowBC4NF
    );
}

// ************************************************************************* //
