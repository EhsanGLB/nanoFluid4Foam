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

#include "thermophysicalFunction4NFFoam.H"
#include "HashTable.H"

#include "colors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermophysicalFunction4NFFoam, 0);
defineRunTimeSelectionTable(thermophysicalFunction4NFFoam, Istream);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<thermophysicalFunction4NFFoam> thermophysicalFunction4NFFoam::New(Istream& is)
{
    if (debug)
    {
        Info<< "thermophysicalFunction4NFFoam::New(Istream&) : "
            << "constructing thermophysicalFunction4NFFoam"
            << endl;
    }

    word thermophysicalFunction4NFFoamType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(thermophysicalFunction4NFFoamType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorIn("thermophysicalFunction4NFFoam::New(Istream&)")
            << "Unknown thermophysicalFunction4NFFoam type "
            << thermophysicalFunction4NFFoamType
            << endl << endl
            << "Valid thermophysicalFunction4NFFoam types are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<thermophysicalFunction4NFFoam>(cstrIter()(is));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
