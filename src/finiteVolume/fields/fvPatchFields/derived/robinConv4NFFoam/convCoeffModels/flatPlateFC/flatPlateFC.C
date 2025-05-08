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

#include "addToRunTimeSelectionTable.H"

#include "dimensionSets4NFFoam.H"//-nanoFluid4Foam
#include "flatPlateFC.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace convCoeffModels
{
    defineTypeNameAndDebug(flatPlateFC, 0);
    addToRunTimeSelectionTable(convCoeffModel, flatPlateFC, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

flatPlateFC::flatPlateFC
(
    const dictionary& dict
)
:
    convCoeffModel(dict),
    Vo_(readScalar(dict.lookup("Vo"))),
    L_(readScalar(dict.lookup("L")))
{
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar flatPlateFC::Re() const
{
    return (Vo_*L_)/nu();
}

scalar flatPlateFC::hoAvg(scalar avgT) const
{
    scalar Re_ = Re();
    scalar hoAvg_ = 0;
    if(Re_ < 5e5)
    {
        hoAvg_ = 0.664*pow(Re_, 0.5)*pow(Pr(), 1/3)*(kappa()/L_);
    }
    else
    {
        hoAvg_ = 0.037*pow(Re_, 0.8)*pow(Pr(), 1/3)*(kappa()/L_);
    }

    return hoAvg_;
}

void flatPlateFC::write(Ostream& os) const
{
    convCoeffModel::write(os);
    os.writeKeyword("Vo") << Vo_ << token::END_STATEMENT << nl;
    os.writeKeyword("L") << L_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace convCoeffModels
} // End namespace Foam

// ************************************************************************* //
