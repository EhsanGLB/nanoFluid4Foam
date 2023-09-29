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

#include "NFDimensionSet.H"//-nanoFluid4Foam
#include "MNEPCM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleModels
{

    defineTypeNameAndDebug(MNEPCM, 0);
    addToRunTimeSelectionTable(particleModel, MNEPCM, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MNEPCM::MNEPCM
(
    const wordList& particleList,
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    particleModel(particleList, nanoFluidPropertiesDict, U, p, T),
    MNEPCMCoreDict_(particleSubDict_.subDict("MNEPCMCoreDict")),
    MNEPCMComponentsc_(MNEPCMCoreDict_.lookup("Components")),
    moleFractionc_(MNEPCMCoreDict_.lookup("moleFraction")),
    materialMixturePtrc_(materialMixture::New(MNEPCMComponentsc_, nanoFluidPropertiesDict)),
    MNEPCMShellDict_(particleSubDict_.subDict("MNEPCMShellDict")),
    MNEPCMComponentss_(MNEPCMShellDict_.lookup("Components")),
    moleFractions_(MNEPCMShellDict_.lookup("moleFraction")),
    materialMixturePtrs_(materialMixture::New(MNEPCMComponentss_, nanoFluidPropertiesDict)),
    weightRatio_("weightRatio", dimless, readScalar(particleSubDict_.lookup("weightRatio"))),
    ds_("d", dimLength, readScalar(particleSubDict_.lookup("d")))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
const volScalarField MNEPCM::rho() const
{
    volScalarField rhoc_(IOobject(("rhoC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoC_"+particleName_), dimDensity, SMALL));
    volScalarField rhos_(IOobject(("rhoS_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoS_"+particleName_), dimDensity, SMALL));
    volScalarField rho_(IOobject(("rho_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rho_"+particleName_), dimDensity, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rhoc_.internalField()[i] = materialMixturePtrc_->rho(pi_, Ti_, moleFractionc_);
        rhos_.internalField()[i] = materialMixturePtrs_->rho(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rhoc_.boundaryField()[i][j] = materialMixturePtrc_->rho(pij_, Tij_, moleFractionc_);
            rhos_.boundaryField()[i][j] = materialMixturePtrs_->rho(pij_, Tij_, moleFractions_);
	}
    }

    rho_ = ( ( 1 + weightRatio_ ) * rhoc_ * rhos_ ) / ( rhos_ + weightRatio_ * rhoc_ );
    
    return rho_;
}


const volScalarField MNEPCM::kappa() const
{
    volScalarField kappac_(IOobject(("kappaC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("kappaC_"+particleName_), dimThermalConductivity, SMALL));
    volScalarField kappas_(IOobject(("kappaS_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("kappaS_"+particleName_), dimThermalConductivity, SMALL));
    volScalarField kappa_(IOobject(("kappa_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("kappa_"+particleName_), dimThermalConductivity, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        kappac_.internalField()[i] = materialMixturePtrc_->K(pi_, Ti_, moleFractionc_);
        kappas_.internalField()[i] = materialMixturePtrs_->K(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            kappac_.boundaryField()[i][j] = materialMixturePtrc_->K(pij_, Tij_, moleFractionc_);
            kappas_.boundaryField()[i][j] = materialMixturePtrs_->K(pij_, Tij_, moleFractions_);
	}
    }

    kappa_ = ds_ / ( ( dc() /kappac_ ) + ( ( ds_ - dc() ) / kappas_ ) );

    return kappa_;
}


const volScalarField MNEPCM::Cp() const
{
    volScalarField Cpc_(IOobject(("CpC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("CpC_"+particleName_), dimSpecificHeatCapacity, SMALL));
    volScalarField hmc_(IOobject(("hmC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("hmC_"+particleName_), dimEnthalpy, SMALL));
    volScalarField Tmrc_(IOobject(("TmrC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("TmrC_"+particleName_), dimTemperature, SMALL));
    volScalarField Cps_(IOobject(("CpS_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("CpS_"+particleName_), dimSpecificHeatCapacity, SMALL));
    volScalarField rhoc_(IOobject(("rhoC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoC_"+particleName_), dimDensity, SMALL));
    volScalarField rhos_(IOobject(("rhoS_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoS_"+particleName_), dimDensity, SMALL));
    volScalarField Cp_(IOobject(("Cp_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Cp_"+particleName_), dimSpecificHeatCapacity, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Cpc_.internalField()[i] = materialMixturePtrc_->cp(pi_, Ti_, moleFractionc_);
        hmc_.internalField()[i] = materialMixturePtrc_->hm(pi_, Ti_, moleFractionc_);
        Tmrc_.internalField()[i] = materialMixturePtrc_->Tmr(pi_, Ti_, moleFractionc_);
        Cps_.internalField()[i] = materialMixturePtrs_->cp(pi_, Ti_, moleFractions_);
        rhoc_.internalField()[i] = materialMixturePtrc_->rho(pi_, Ti_, moleFractionc_);
        rhos_.internalField()[i] = materialMixturePtrs_->rho(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Cpc_.boundaryField()[i][j] = materialMixturePtrc_->cp(pij_, Tij_, moleFractionc_);
            hmc_.boundaryField()[i][j] = materialMixturePtrc_->hm(pij_, Tij_, moleFractionc_);
            Tmrc_.boundaryField()[i][j] = materialMixturePtrc_->Tmr(pij_, Tij_, moleFractionc_);
            Cps_.boundaryField()[i][j] = materialMixturePtrs_->cp(pij_, Tij_, moleFractions_);
            rhoc_.boundaryField()[i][j] = materialMixturePtrc_->rho(pij_, Tij_, moleFractionc_);
            rhos_.boundaryField()[i][j] = materialMixturePtrs_->rho(pij_, Tij_, moleFractions_);
	}
    }

    volScalarField CpC_ = Cpc_ + (3.14/2) * ( ( hmc_ / Tmrc_ ) - Cpc_ ) * fc();
    Cp_ = ( ( CpC_ + weightRatio_ * Cps_) * rhoc_ * rhos_ ) / ( ( rhos_ + weightRatio_ * rhoc_ ) * rho() );

    return Cp_;
}


const volScalarField MNEPCM::beta() const
{


    volScalarField rhoc_(IOobject(("rhoC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoC_"+particleName_), dimDensity, SMALL));
    volScalarField rhos_(IOobject(("rhoS_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoS_"+particleName_), dimDensity, SMALL));
    volScalarField betac_(IOobject(("betaC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("betaC_"+particleName_), dimThermalExpansion, SMALL));
    volScalarField betas_(IOobject(("betaS_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("betaS_"+particleName_), dimThermalExpansion, SMALL));
    volScalarField beta_(IOobject(("beta_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("beta_"+particleName_), dimThermalExpansion, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rhoc_.internalField()[i] = materialMixturePtrc_->rho(pi_, Ti_, moleFractionc_);
        rhos_.internalField()[i] = materialMixturePtrs_->rho(pi_, Ti_, moleFractions_);
        betac_.internalField()[i] = materialMixturePtrc_->beta(pi_, Ti_, moleFractionc_);
        betas_.internalField()[i] = materialMixturePtrs_->beta(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rhoc_.boundaryField()[i][j] = materialMixturePtrc_->rho(pij_, Tij_, moleFractionc_);
            rhos_.boundaryField()[i][j] = materialMixturePtrs_->rho(pij_, Tij_, moleFractions_);
            betac_.boundaryField()[i][j] = materialMixturePtrc_->beta(pij_, Tij_, moleFractionc_);
            betas_.boundaryField()[i][j] = materialMixturePtrs_->beta(pij_, Tij_, moleFractions_);
	}
    }

    beta_ = betac_ + ( ( betas_ - betac_ ) / 2 ) * ( 1 - weightRatio_ * rhos_ / rhoc_ );

    return beta_;
}


const volScalarField MNEPCM::rhoR() const
{
    volScalarField rhoRc_(IOobject(("rhoRc_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoRc_"+particleName_), dimElectricalResistivity, SMALL));
    volScalarField rhoRs_(IOobject(("rhoRs_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoRs_"+particleName_), dimElectricalResistivity, SMALL));
    volScalarField rhoR_(IOobject(("rhoR_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoR_"+particleName_), dimElectricalResistivity, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rhoRc_.internalField()[i] = materialMixturePtrc_->rhoR(pi_, Ti_, moleFractionc_);
        rhoRs_.internalField()[i] = materialMixturePtrs_->rhoR(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rhoRc_.boundaryField()[i][j] = materialMixturePtrc_->rhoR(pij_, Tij_, moleFractionc_);
            rhoRs_.boundaryField()[i][j] = materialMixturePtrs_->rhoR(pij_, Tij_, moleFractions_);
	}
    }

    rhoR_ = ( ( dc() * rhoRc_ ) + ( ( ds_ - dc() ) * rhoRs_ ) ) / ds_;

    return rhoR_;
}


const volScalarField MNEPCM::mu() const
{
    volScalarField muc_(IOobject(("muc_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("muc_"+particleName_), dimDynamicViscosity, SMALL));
    volScalarField mus_(IOobject(("mus_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("mus_"+particleName_), dimDynamicViscosity, SMALL));
    volScalarField mu_(IOobject(("mu_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("mu_"+particleName_), dimDynamicViscosity, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        muc_.internalField()[i] = materialMixturePtrc_->mu(pi_, Ti_, moleFractionc_);
        mus_.internalField()[i] = materialMixturePtrs_->mu(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            muc_.boundaryField()[i][j] = materialMixturePtrc_->mu(pij_, Tij_, moleFractionc_);
            mus_.boundaryField()[i][j] = materialMixturePtrs_->mu(pij_, Tij_, moleFractions_);
	}
    }
    
    mu_ = mus_;

    return mu_;
}


const volScalarField MNEPCM::pv() const
{
    volScalarField pvc_(IOobject(("pvc_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("pvc_"+particleName_), dimPressure, SMALL));
    volScalarField pvs_(IOobject(("pvs_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("pvs_"+particleName_), dimPressure, SMALL));
    volScalarField pv_(IOobject(("pv_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("pv_"+particleName_), dimPressure, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        pvc_.internalField()[i] = materialMixturePtrc_->pv(pi_, Ti_, moleFractionc_);
        pvs_.internalField()[i] = materialMixturePtrs_->pv(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            pvc_.boundaryField()[i][j] = materialMixturePtrc_->pv(pij_, Tij_, moleFractionc_);
            pvs_.boundaryField()[i][j] = materialMixturePtrs_->pv(pij_, Tij_, moleFractions_);
	}
    }
    
    pv_ = pvs_;

    return pv_;
}


const volScalarField MNEPCM::hm() const
{
    volScalarField hmc_(IOobject(("hmc_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("hmc_"+particleName_), dimEnthalpy, SMALL));
    volScalarField hms_(IOobject(("hms_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("hms_"+particleName_), dimEnthalpy, SMALL));
    volScalarField hm_(IOobject(("hm_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("hm_"+particleName_), dimEnthalpy, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        hmc_.internalField()[i] = materialMixturePtrc_->hm(pi_, Ti_, moleFractionc_);
        hms_.internalField()[i] = materialMixturePtrs_->hm(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            hmc_.boundaryField()[i][j] = materialMixturePtrc_->hm(pij_, Tij_, moleFractionc_);
            hms_.boundaryField()[i][j] = materialMixturePtrs_->hm(pij_, Tij_, moleFractions_);
	}
    }
    
    hm_ = hms_;

    return hm_;
}


const volScalarField MNEPCM::Tm() const
{
    volScalarField Tmc_(IOobject(("Tmc_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Tmc_"+particleName_), dimTemperature, SMALL));
    volScalarField Tms_(IOobject(("Tms_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Tms_"+particleName_), dimTemperature, SMALL));
    volScalarField Tm_(IOobject(("Tm_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Tm_"+particleName_), dimTemperature, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Tmc_.internalField()[i] = materialMixturePtrc_->Tm(pi_, Ti_, moleFractionc_);
        Tms_.internalField()[i] = materialMixturePtrs_->Tm(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Tmc_.boundaryField()[i][j] = materialMixturePtrc_->Tm(pij_, Tij_, moleFractionc_);
            Tms_.boundaryField()[i][j] = materialMixturePtrs_->Tm(pij_, Tij_, moleFractions_);
	}
    }
    
    Tm_ = Tms_;

    return Tm_;
}


const volScalarField MNEPCM::Tmr() const
{
    volScalarField Tmrc_(IOobject(("Tmrc_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Tmrc_"+particleName_), dimTemperature, SMALL));
    volScalarField Tmrs_(IOobject(("Tmrs_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Tmrs_"+particleName_), dimTemperature, SMALL));
    volScalarField Tmr_(IOobject(("Tmr_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Tmr_"+particleName_), dimTemperature, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Tmrc_.internalField()[i] = materialMixturePtrc_->Tmr(pi_, Ti_, moleFractionc_);
        Tmrs_.internalField()[i] = materialMixturePtrs_->Tmr(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Tmrc_.boundaryField()[i][j] = materialMixturePtrc_->Tmr(pij_, Tij_, moleFractionc_);
            Tmrs_.boundaryField()[i][j] = materialMixturePtrs_->Tmr(pij_, Tij_, moleFractions_);
	}
    }
    
    Tmr_ = Tmrs_;

    return Tmr_;
}


const volScalarField MNEPCM::molten() const
{
    volScalarField Tmc_(IOobject(("TmC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("TmC_"+particleName_), dimTemperature, SMALL));
    volScalarField Tmrc_(IOobject(("TmrC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("TmrC_"+particleName_), dimTemperature, SMALL));
    volScalarField molten_(IOobject("molten", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("molten", dimless, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Tmc_.internalField()[i] = materialMixturePtrc_->Tm(pi_, Ti_, moleFractionc_);
        Tmrc_.internalField()[i] = materialMixturePtrc_->Tmr(pi_, Ti_, moleFractionc_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Tmc_.boundaryField()[i][j] = materialMixturePtrc_->Tm(pij_, Tij_, moleFractionc_);
            Tmrc_.boundaryField()[i][j] = materialMixturePtrc_->Tmr(pij_, Tij_, moleFractionc_);
	}
    }

	forAll(T_, celli)
	{
		if (T_[celli] < scalar( Tmc_[celli] - Tmrc_[celli]/2 ) )
		{
			molten_[celli] = scalar(0.0);
		}
		else if ( scalar( Tmc_[celli] - Tmrc_[celli]/2 ) <= T_[celli] && T_[celli] <= scalar( Tmc_[celli] + Tmrc_[celli]/2) )
		{
			molten_[celli] = sin ( ( 3.14 / Tmrc_[celli] ) * ( T_[celli] - ( Tmc_[celli] - Tmrc_[celli] / 2 ) ) );
		}
		else if (T_[celli] > scalar(Tmc_[celli] + Tmrc_[celli]/2 ) )
		{
			molten_[celli] = scalar(1.0);
		}
	}

    return molten_;
}


const volScalarField MNEPCM::fc() const
{
    volScalarField Tmc_(IOobject(("TmC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("TmC_"+particleName_), dimTemperature, SMALL));
    volScalarField Tmrc_(IOobject(("TmrC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("TmrC_"+particleName_), dimTemperature, SMALL));
    volScalarField fc_(IOobject(("fc_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("fc_"+particleName_), dimless, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Tmc_.internalField()[i] = materialMixturePtrc_->Tm(pi_, Ti_, moleFractionc_);
        Tmrc_.internalField()[i] = materialMixturePtrc_->Tmr(pi_, Ti_, moleFractionc_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Tmc_.boundaryField()[i][j] = materialMixturePtrc_->Tm(pij_, Tij_, moleFractionc_);
            Tmrc_.boundaryField()[i][j] = materialMixturePtrc_->Tmr(pij_, Tij_, moleFractionc_);
	}
    }


	forAll(T_, celli)
	{
		if (T_[celli] < scalar( Tmc_[celli] - Tmrc_[celli]/2 ) )
		{
			fc_[celli] = scalar(0.0);
		}
		else if ( scalar( Tmc_[celli] - Tmrc_[celli]/2 ) <= T_[celli] && T_[celli] <= scalar( Tmc_[celli] + Tmrc_[celli]/2) )
		{
			fc_[celli] = sin ( ( 3.14 / Tmrc_[celli] ) * ( T_[celli] - ( Tmc_[celli] - Tmrc_[celli] / 2 ) ) );
		}
		else if (T_[celli] > scalar(Tmc_[celli] + Tmrc_[celli]/2 ) )
		{
			fc_[celli] = scalar(0.0);
		}
	}

    return fc_;
}


const volScalarField MNEPCM::dc() const
{
    volScalarField rhoc_(IOobject(("rhoC_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoC_"+particleName_), dimDensity, SMALL));
    volScalarField rhos_(IOobject(("rhoS_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoS_"+particleName_), dimDensity, SMALL));
    volScalarField dc_(IOobject(("dc_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("dc_"+particleName_), dimLength, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rhoc_.internalField()[i] = materialMixturePtrc_->rho(pi_, Ti_, moleFractionc_);
        rhos_.internalField()[i] = materialMixturePtrs_->rho(pi_, Ti_, moleFractions_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rhoc_.boundaryField()[i][j] = materialMixturePtrc_->rho(pij_, Tij_, moleFractionc_);
            rhos_.boundaryField()[i][j] = materialMixturePtrs_->rho(pij_, Tij_, moleFractions_);
	}
    }

    dc_ = ds_ * pow ( rhos_ / ( rhos_ + weightRatio_ * rhoc_ ) , scalar( 0.3333 ) );

    return dc_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace particleModels
} // End namespace Foam

// ************************************************************************* //
