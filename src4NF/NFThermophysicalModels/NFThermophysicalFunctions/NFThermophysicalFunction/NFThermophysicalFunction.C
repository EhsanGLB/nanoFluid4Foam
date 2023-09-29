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

#include "error.H"

#include "NFThermophysicalFunction.H"
#include "HashTable.H"

#include "colors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(NFThermophysicalFunction, 0);
defineRunTimeSelectionTable(NFThermophysicalFunction, Istream);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<NFThermophysicalFunction> NFThermophysicalFunction::New(Istream& is)
{
    if (debug)
    {
        Info<< "NFThermophysicalFunction::New(Istream&) : "
            << "constructing NFThermophysicalFunction"
            << endl;
    }

    word NFThermophysicalFunctionType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(NFThermophysicalFunctionType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorIn("NFThermophysicalFunction::New(Istream&)")
            << "Unknown NFThermophysicalFunction type "
            << NFThermophysicalFunctionType
            << endl << endl
            << "Valid NFThermophysicalFunction types are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<NFThermophysicalFunction>(cstrIter()(is));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
