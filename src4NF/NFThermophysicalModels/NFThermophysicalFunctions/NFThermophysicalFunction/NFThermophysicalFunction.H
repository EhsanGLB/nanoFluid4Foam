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

Class
    Foam::NFThermophysicalFunction

Description
    Abstract base class for thermo-physical functions

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef NFThermophysicalFunction_H
#define NFThermophysicalFunction_H

#include "scalar.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class NFThermophysicalFunction Declaration
\*---------------------------------------------------------------------------*/

class NFThermophysicalFunction
{

public:

    //- Runtime type information
    TypeName("NFThermophysicalFunction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            NFThermophysicalFunction,
            Istream,
            (Istream& is),
            (is)
        );


    // Constructors

        //- Construct null
        NFThermophysicalFunction()
        {}

        //- Return a pointer to a new NFThermophysicalFunction created from input
        static autoPtr<NFThermophysicalFunction> New(Istream& is);


    //- Destructor
    virtual ~NFThermophysicalFunction()
    {}


    // Member Functions

        //- Evaluate the function and return the result
        virtual scalar f(scalar T) const = 0;

        //- Write the function coefficients
        virtual void writeData(Ostream& os) const = 0;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const NFThermophysicalFunction& f)
        {
            f.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
