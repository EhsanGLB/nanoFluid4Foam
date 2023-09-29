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
    Foam::complex

Description
    Extension to the c++ complex library type.

SourceFiles
    complexI.H
    complex.C

\*---------------------------------------------------------------------------*/

//#ifndef complex_H
//#define complex_H

#include "scalar.H"
#include "bool.H"
#include "word.H"
#include "contiguous.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline unsigned long int factorial(unsigned long int n);

// Forward declaration of friend functions and operators
inline complex exp(const complex&);
inline complex log(const complex&);
inline complex pow(const complex&, const scalar);
inline complex sqrt(const complex&);
inline complex sin(const complex&);
inline complex cos(const complex&);
inline complex J(const int&, const complex&);
/*---------------------------------------------------------------------------*\
                           Class complex Declaration
\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "complexI4NF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#endif

// ************************************************************************* //
