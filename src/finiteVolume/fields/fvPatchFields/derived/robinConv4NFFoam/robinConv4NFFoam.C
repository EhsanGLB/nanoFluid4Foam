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

#include "robinConv4NFFoam.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::robinConv4NFFoam::robinConv4NFFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    dictc_(),
    TName_("T"),
    kappaEffName_("kappaEff"),
    fluxLoc_(0.0),
    fluxTime_(0.0),
    dir_(0, 0, 0),
    To_(0.0),
    convCoeffModelPtr_()
{}


Foam::robinConv4NFFoam::robinConv4NFFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    dictc_(dict),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaEffName_(dict.lookupOrDefault<word>("kappaEff", "kappaEff")),
    fluxLoc_(),
    fluxTime_(),
    dir_(dict.lookup("dir")),
    To_(readScalar(dict.lookup("To"))),
    convCoeffModelPtr_(convCoeffModel::New(dict))
{
    Istream& isfluxLoc_ = dict.lookup("fluxLoc");
    isfluxLoc_.format(IOstream::ASCII);
    isfluxLoc_ >> fluxLoc_;

    Istream& isfluxTime_ = dict.lookup("fluxTime");
    isfluxTime_.format(IOstream::ASCII);
    isfluxTime_ >> fluxTime_;

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


Foam::robinConv4NFFoam::robinConv4NFFoam
(
    const robinConv4NFFoam& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    dictc_(ptf.dictc_),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    fluxLoc_(ptf.fluxLoc_),
    fluxTime_(ptf.fluxTime_),
    dir_(ptf.dir_),
    To_(ptf.To_),
    convCoeffModelPtr_(ptf.convCoeffModelPtr_)
{}


Foam::robinConv4NFFoam::robinConv4NFFoam
(
    const robinConv4NFFoam& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    dictc_(ptf.dictc_),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    fluxLoc_(ptf.fluxLoc_),
    fluxTime_(ptf.fluxTime_),
    dir_(ptf.dir_),
    To_(ptf.To_),
    convCoeffModelPtr_(ptf.convCoeffModelPtr_)
{}


Foam::robinConv4NFFoam::robinConv4NFFoam
(
    const robinConv4NFFoam& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    dictc_(ptf.dictc_),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    fluxLoc_(ptf.fluxLoc_),
    fluxTime_(ptf.fluxTime_),
    dir_(ptf.dir_),
    To_(ptf.To_),
    convCoeffModelPtr_(ptf.convCoeffModelPtr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::robinConv4NFFoam::hoAvg()
{
    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    scalar avgT_ = sum(Tp*patch().magSf())/sum(patch().magSf());
    return convCoeffModelPtr_->hoAvg(avgT_);
}

Foam::scalar Foam::robinConv4NFFoam::fluxLocFunction(scalar loc)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < fluxLoc_.size() )
    {
        f_ += fluxLoc_[n_]*pow(loc, n_);
        n_ += 1;
    }

    return f_;
}


Foam::scalar Foam::robinConv4NFFoam::fluxTimeFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < fluxTime_.size() )
    {
        f_ += fluxTime_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}

void Foam::robinConv4NFFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- get the currernt time and time step
    scalar t_ = this->db().time().value();


    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& kappaEffp = lookupPatchField<volScalarField, scalar>(kappaEffName_);
    scalarField patchLoc_ = patch().Cf() & dir_;


    //scalar len_ = max(patchLoc_) - min(patchLoc_);
    Info << "hoAvg_: " << hoAvg() << endl;


    forAll(patchLoc_, i)
    {
        gradient()[i] = ( fluxTimeFunction(t_)*fluxLocFunction(patchLoc_[i]) - hoAvg() * (Tp[i] - To_) )/kappaEffp[i];
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::robinConv4NFFoam::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    convCoeffModelPtr_->write(os);
    os.writeKeyword("fluxLoc") << fluxLoc_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluxTime") << fluxTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("dir") << dir_ << token::END_STATEMENT << nl;
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
        robinConv4NFFoam
    );
}

// ************************************************************************* //
