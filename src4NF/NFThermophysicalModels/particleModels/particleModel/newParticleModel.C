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
#include "surfaceFields.H"

#include "particleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<particleModel> particleModel::New
(
    const wordList& particleList,
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
{
    word particleModelTypeName(particleList[1]);

    Info<< "Selecting particle model "
        << particleModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(particleModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "particleModel::New(const wordList&, "
	    "const dictionary&, "
	    "const volVectorField&, "
	    "const volScalarField&, "
            "const volScalarField&)"
        )   << "Unknown particleModel type "
            << particleModelTypeName << endl << endl
            << "Valid  particleModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<particleModel>
        (cstrIter()(particleList, nanoFluidPropertiesDict, U, p, T));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
