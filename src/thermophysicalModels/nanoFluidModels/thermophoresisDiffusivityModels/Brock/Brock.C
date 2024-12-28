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
#include "Brock.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophoresisDiffusivityModels
{
    defineTypeNameAndDebug(Brock, 0);
    addToRunTimeSelectionTable(thermophoresisDiffusivityModel, Brock, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Brock::Brock
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
    thermophoresisDiffusivityModel(nanoFluidPropertiesDict, U, p, T, alphasPtr, baseFluidPtr, particlesProperties),
    BrockCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    Cs_(readScalar(BrockCoeffs_.lookup("Cs"))),
    Cm_(readScalar(BrockCoeffs_.lookup("Cm"))),
    Ct_(readScalar(BrockCoeffs_.lookup("Ct"))),
    lbf_(readScalar(BrockCoeffs_.lookup("lbf"))),
    dnp_(BrockCoeffs_.lookup("dnp")),
    DTs_(particlesProperties.size())
{
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField mubf_ = baseFluidPtr_->mu();

    forAll(DTs_, i)
    {
        DTs_.set
        (
            i,
            new volScalarField (IOobject(("DT_"+particlesProperties_[i].particleName()), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("DT_"+particlesProperties_[i].particleName()), dimDiffusivity, SMALL))
        );
    }

    forAll(DTs_, i)
    {
        scalar Kn_ = lbf_ / dnp_[i] ;
        volScalarField beta_ = ( 2.0 * Cs_ * ( kappabf_ + particlesProperties_[i].kappa() * Kn_ ) ) / ( ( 1.0 + 3.0 * Cm_ * Kn_ ) * ( 2.0 * kappabf_ + particlesProperties_[i].kappa() + 2.0 * particlesProperties_[i].kappa() * Ct_ * Kn_ ) );
        DTs_[i] = beta_ * ( mubf_ / rhobf_ ) * alphasPtr_[i];
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermophoresisDiffusivityModels
} // End namespace Foam

// ************************************************************************* //
