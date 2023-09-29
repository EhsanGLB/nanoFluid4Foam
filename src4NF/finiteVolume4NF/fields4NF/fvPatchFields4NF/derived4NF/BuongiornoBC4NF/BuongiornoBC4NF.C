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

#include "BuongiornoBC4NF.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BuongiornoBC4NF::BuongiornoBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    alphaName_("alpha"),
    particleName_("alpha"),
    massFlux_(0.0),
    ho_(0.0),
    alphao_(0.0)
{}


Foam::BuongiornoBC4NF::BuongiornoBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    particleName_(dict.lookup("name")),
    massFlux_(readScalar(dict.lookup("massFlux"))),
    ho_(readScalar(dict.lookup("ho"))),
    alphao_(readScalar(dict.lookup("alphao")))
{
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


Foam::BuongiornoBC4NF::BuongiornoBC4NF
(
    const BuongiornoBC4NF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    particleName_(ptf.particleName_),
    massFlux_(ptf.massFlux_),
    ho_(ptf.ho_),
    alphao_(ptf.alphao_)
{}


Foam::BuongiornoBC4NF::BuongiornoBC4NF
(
    const BuongiornoBC4NF& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    particleName_(ptf.particleName_),
    massFlux_(ptf.massFlux_),
    ho_(ptf.ho_),
    alphao_(ptf.alphao_)
{}


Foam::BuongiornoBC4NF::BuongiornoBC4NF
(
    const BuongiornoBC4NF& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    particleName_(ptf.particleName_),
    massFlux_(ptf.massFlux_),
    ho_(ptf.ho_),
    alphao_(ptf.alphao_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BuongiornoBC4NF::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& alphap = lookupPatchField<volScalarField, scalar>("alpha_" + particleName_);
    const fvPatchScalarField& DBEffp = lookupPatchField<volScalarField, scalar>("DBEff_" + particleName_);
    const fvPatchScalarField& DTp = lookupPatchField<volScalarField, scalar>("DT_" + particleName_);

    gradient() = massFlux_ - ho_ * (alphap - alphao_) - ((DTp/Tp)/DBEffp)*Tp.snGrad();
    fixedGradientFvPatchScalarField::updateCoeffs();

}


void Foam::BuongiornoBC4NF::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("name") << particleName_ << token::END_STATEMENT << nl;
    os.writeKeyword("massFlux") << massFlux_ << token::END_STATEMENT << nl;
    os.writeKeyword("ho") << ho_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphao") << alphao_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        BuongiornoBC4NF
    );
}

// ************************************************************************* //
