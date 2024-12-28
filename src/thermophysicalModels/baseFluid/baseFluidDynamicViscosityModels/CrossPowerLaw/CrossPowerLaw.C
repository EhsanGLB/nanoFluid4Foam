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

#include "dimensionSets4NFFoam.H"
#include "CrossPowerLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace baseFluidDynamicViscosityModels
{
    defineTypeNameAndDebug(CrossPowerLaw, 0);
    addToRunTimeSelectionTable(baseFluidDynamicViscosityModel, CrossPowerLaw, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CrossPowerLaw::CrossPowerLaw
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    baseFluidDynamicViscosityModel(nanoFluidPropertiesDict, U, p, T),
    CrossPowerLawCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    mu0_(CrossPowerLawCoeffs_.lookup("mu0")),
    muInf_(CrossPowerLawCoeffs_.lookup("muInf")),
    m_(CrossPowerLawCoeffs_.lookup("m")),
    n_(CrossPowerLawCoeffs_.lookup("n"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField CrossPowerLaw::mu() const
{
    volScalarField mu_(IOobject("muBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muBF", dimDynamicViscosity, SMALL));

    mu_ = (mu0_ - muInf_)/(scalar(1) + pow(m_*strainRate(), n_)) + muInf_;

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace baseFluidDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
