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
#include "fvcGrad.H"

#include "baseFluidDynamicViscosityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(baseFluidDynamicViscosityModel, 0);
    defineRunTimeSelectionTable(baseFluidDynamicViscosityModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::baseFluidDynamicViscosityModel::baseFluidDynamicViscosityModel
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    nanoFluidPropertiesDict_(nanoFluidPropertiesDict),
    U_(U),
    p_(p),
    T_(T)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::baseFluidDynamicViscosityModel::strainRate() const
{
    // Bug fix: sqrt(2) inconsistency.  HJ, 8/Dec/2009
    return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}

} // End namespace Foam

// ************************************************************************* //
