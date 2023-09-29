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

#include "RobinBC4NF.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RobinBC4NF::RobinBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    kappaEffName_("kappaEff"),
    flux_(0.0),
    dir_(0, 0, 0),
    ho_(0.0),
    To_(0.0)
{}


Foam::RobinBC4NF::RobinBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaEffName_(dict.lookupOrDefault<word>("kappaEff", "kappaEff")),
    flux_(),
    dir_(dict.lookup("dir")),
    ho_(readScalar(dict.lookup("ho"))),
    To_(readScalar(dict.lookup("To")))
{
    Istream& isflux_ = dict.lookup("flux");
    isflux_.format(IOstream::ASCII);
    isflux_ >> flux_;

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


Foam::RobinBC4NF::RobinBC4NF
(
    const RobinBC4NF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    flux_(ptf.flux_),
    dir_(ptf.dir_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


Foam::RobinBC4NF::RobinBC4NF
(
    const RobinBC4NF& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    flux_(ptf.flux_),
    dir_(ptf.dir_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


Foam::RobinBC4NF::RobinBC4NF
(
    const RobinBC4NF& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    flux_(ptf.flux_),
    dir_(ptf.dir_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::RobinBC4NF::fluxFunction(scalar loc)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < flux_.size() )
    {
        f_ += flux_[n_]*pow(loc, n_);
        n_ += 1;
    }

    return f_;
}


void Foam::RobinBC4NF::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& kappaEffp = lookupPatchField<volScalarField, scalar>(kappaEffName_);
    scalarField patchLoc_ = patch().Cf() & dir_;

    forAll(patchLoc_, i)
    {
        gradient()[i] = ( fluxFunction(patchLoc_[i]) - ho_ * (Tp[i] - To_) )/kappaEffp[i];
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::RobinBC4NF::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("flux") << flux_ << token::END_STATEMENT << nl;
    os.writeKeyword("dir") << dir_ << token::END_STATEMENT << nl;
    os.writeKeyword("ho") << ho_ << token::END_STATEMENT << nl;
    os.writeKeyword("To") << To_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        RobinBC4NF
    );
}

// ************************************************************************* //
