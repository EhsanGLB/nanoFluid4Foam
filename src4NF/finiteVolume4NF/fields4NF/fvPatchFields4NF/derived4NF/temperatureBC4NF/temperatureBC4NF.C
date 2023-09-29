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

#include "temperatureBC4NF.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureBC4NF::temperatureBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    temperature_(0.0),
    dir_(0, 0, 0)
{}


Foam::temperatureBC4NF::temperatureBC4NF
(
    const temperatureBC4NF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    temperature_(ptf.temperature_),
    dir_(ptf.dir_)
{}


Foam::temperatureBC4NF::temperatureBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    temperature_(),
    dir_(dict.lookup("dir"))
{
    Istream& istemperature_ = dict.lookup("temperature");
    istemperature_.format(IOstream::ASCII);
    istemperature_ >> temperature_;
}


Foam::temperatureBC4NF::temperatureBC4NF
(
    const temperatureBC4NF& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    temperature_(fcvpvf.temperature_),
    dir_(fcvpvf.dir_)
{}


Foam::temperatureBC4NF::temperatureBC4NF
(
    const temperatureBC4NF& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    temperature_(ptf.temperature_),
    dir_(ptf.dir_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::temperatureBC4NF::temperatureFunction(scalar loc)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < temperature_.size() )
    {
        f_ += temperature_[n_]*pow(loc, n_);
        n_ += 1;
    }

    return f_;
}


void Foam::temperatureBC4NF::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    // Get range and orientation
    boundBox bb(patch().patch().localPoints(), true);
    vector ctr = 0.5*(bb.max() + bb.min()) - 0.5*(bb.max() - bb.min());
    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate
    scalarField coord = (c - ctr) & dir_;
    scalarField valueTemp(coord.size());

    forAll(valueTemp, i)
    {
        valueTemp[i] = temperatureFunction(coord[i]);
    }

    scalarField::operator=(valueTemp);
    fixedValueFvPatchScalarField::updateCoeffs();
}


// Write
void Foam::temperatureBC4NF::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("temperature") << temperature_ << token::END_STATEMENT << nl;
    os.writeKeyword("dir") << dir_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField(fvPatchScalarField, temperatureBC4NF);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// ************************************************************************* //
