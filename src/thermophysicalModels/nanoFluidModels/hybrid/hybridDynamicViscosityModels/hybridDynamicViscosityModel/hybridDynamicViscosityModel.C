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

#include "hybridDynamicViscosityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hybridDynamicViscosityModel, 0);
    defineRunTimeSelectionTable(hybridDynamicViscosityModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hybridDynamicViscosityModel::hybridDynamicViscosityModel
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T,
    const PtrList<volScalarField>& alphasPtr,
    const autoPtr<baseFluid>& baseFluidPtr,
    const PtrList<particleModel>& particlesProperties
)
:
    nanoFluidPropertiesDict_(nanoFluidPropertiesDict),
    U_(U),
    p_(p),
    T_(T),
    alphasPtr_(alphasPtr),
    baseFluidPtr_(baseFluidPtr),
    particlesProperties_(particlesProperties),
    particlesComponents_(particlesProperties_.size()),
    kB_("kB", dimensionSet(1, 2, -2, -1, 0, 0, 0), 1.38064e-23)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
