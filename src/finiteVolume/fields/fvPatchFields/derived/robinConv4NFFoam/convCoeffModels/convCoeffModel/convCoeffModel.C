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

#include "volFields.H"

#include "dimensionSets4NFFoam.H"//-nanoFluid4Foam
#include "convCoeffModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(convCoeffModel, 0);
    defineRunTimeSelectionTable(convCoeffModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

convCoeffModel::convCoeffModel
(
    const dictionary& dict
)
:
    dict_(dict),
    convCoeffModelType_(dict.lookup("convCoeffModel")),
    To_(readScalar(dict.lookup("To"))),
    flowProp_()
{
    Istream& isflowProp_ = dict.lookup("flowProp");
    isflowProp_.format(IOstream::ASCII);
    isflowProp_ >> flowProp_;
    rhoc_ = flowProp_[0];
    kappac_ = flowProp_[1];
    Cpc_ = flowProp_[2];
    muc_ = flowProp_[3];
    betac_ = flowProp_[4];
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar convCoeffModel::rho() const
{
    return rhoc_;
}

scalar convCoeffModel::kappa() const
{
    return kappac_;
}

scalar convCoeffModel::Cp() const
{
    return Cpc_;
}

scalar convCoeffModel::mu() const
{
    return muc_;
}

scalar convCoeffModel::beta() const
{
    return betac_;
}

scalar convCoeffModel::nu() const
{
    return muc_/rhoc_;
}

scalar convCoeffModel::alpha() const
{
    return kappac_/(rhoc_*Cpc_);
}

scalar convCoeffModel::Pr() const
{
    return nu()/alpha();
}

void convCoeffModel::write(Ostream& os) const
{
    os.writeKeyword("convCoeffModel") << convCoeffModelType_ << token::END_STATEMENT << nl;
    os.writeKeyword("flowProp") << flowProp_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
