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

#include "Devi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace hybridThermalConductivityModels
{
    defineTypeNameAndDebug(Devi, 0);
    addToRunTimeSelectionTable(hybridThermalConductivityModel, Devi, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Devi::Devi
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
    hybridThermalConductivityModel(nanoFluidPropertiesDict, U, p, T, alphasPtr, baseFluidPtr, particlesProperties),
    DeviCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    lambdaList_(DeviCoeffs_.lookup("lambda"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Devi::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();

    kappa_ = kappabf_;
    forAll(particlesComponents_, i)
    {
        kappa_ *= ( ( particlesProperties_[i].kappa() + ( lambdaList_[i] - 1 ) * kappabf_ - ( lambdaList_[i] - 1 ) * ( kappabf_ - particlesProperties_[i].kappa() ) * alphasPtr_[i] ) / ( particlesProperties_[i].kappa() + ( lambdaList_[i] - 1 ) * kappabf_ + ( kappabf_ - particlesProperties_[i].kappa() ) * alphasPtr_[i] ) );
    }

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hybridThermalConductivityModel
} // End namespace Foam

// ************************************************************************* //
