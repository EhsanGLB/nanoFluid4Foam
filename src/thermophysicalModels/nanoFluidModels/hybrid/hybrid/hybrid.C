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
#include "hybrid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nanoFluidModels
{
    defineTypeNameAndDebug(hybrid, 0);
    addToRunTimeSelectionTable(nanoFluidModel, hybrid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hybrid::hybrid
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T,
    const PtrList<volScalarField>& alphasPtr
)
:
    nanoFluidModel(nanoFluidPropertiesDict, U, p, T, alphasPtr),
    hybridThermalConductivityModelPtr_(hybridThermalConductivityModel::New(nanoFluidPropertiesDict, U, p, T, alphasPtr_, baseFluidPtr_, particlesProperties_)),
    hybridDynamicViscosityModelPtr_(hybridDynamicViscosityModel::New(nanoFluidPropertiesDict, U, p, T, alphasPtr_, baseFluidPtr_, particlesProperties_))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField hybrid::rho() const
{
    volScalarField rho_(IOobject("rhoNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("rhoNF", dimDensity, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();

    rho_ = ( 1 - alpha() ) * rhobf_;
    forAll(particlesComponents_, i)
    {
        rho_ += alphasPtr_[i] * particlesProperties_[i].rho();
    }

    return rho_;
}


const volScalarField hybrid::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();

    kappa_ = hybridThermalConductivityModelPtr_->kappa();;

    return kappa_;
}


const volScalarField hybrid::Cp() const
{
    volScalarField Cp_(IOobject("CpNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("CpNF", dimSpecificHeatCapacity, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField Cpbf_ = baseFluidPtr_->Cp();

    Cp_ = ( ( 1 - alpha() ) * rhobf_ * Cpbf_ ) / rho();
    forAll(particlesComponents_, i)
    {
        Cp_ += ( alphasPtr_[i] * particlesProperties_[i].rho() * particlesProperties_[i].Cp() ) / rho();
    }

    return Cp_;
}


const volScalarField hybrid::beta() const
{
    volScalarField beta_(IOobject("betaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("betaNF", dimThermalExpansion, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField betabf_ = baseFluidPtr_->beta();

    beta_ = ( ( 1 - alpha() ) * rhobf_ * betabf_ ) / rho();
    forAll(particlesComponents_, i)
    {
        beta_ += ( alphasPtr_[i] * particlesProperties_[i].rho() * particlesProperties_[i].beta() ) / rho();
    }

    return beta_;
}


const volScalarField hybrid::rhoR() const
{
    volScalarField rhoR_(IOobject("rhoRNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("rhoRNF", dimElectricalResistivity, SMALL));
    volScalarField rhoRbf_ = baseFluidPtr_->rhoR();

    rhoR_ = rhoRbf_;

    return rhoR_;
}


const volScalarField hybrid::mu() const
{
    volScalarField mu_(IOobject("muNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muNF", dimDynamicViscosity, SMALL));
    volScalarField mubf_ = baseFluidPtr_->mu();

    mu_ = hybridDynamicViscosityModelPtr_->mu();

    return mu_;
}


const volScalarField hybrid::pv() const
{
    volScalarField pv_(IOobject("pvNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("pvNF", dimPressure, SMALL));
    volScalarField pvbf_ = baseFluidPtr_->pv();

    pv_ = pvbf_;

    return pv_;
}


const volScalarField hybrid::hm() const
{
    volScalarField hm_(IOobject("hmNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("hmNF", dimEnthalpy, SMALL));
    volScalarField hmbf_ = baseFluidPtr_->hm();

    hm_ = ( 1 - alpha() ) * hmbf_;

    return hm_;
}


const volScalarField hybrid::Tm() const
{
    volScalarField Tm_(IOobject("TmNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("TmNF", dimTemperature, SMALL));
    volScalarField Tmbf_ = baseFluidPtr_->Tm();

    Tm_ = Tmbf_;

    return Tm_;
}


const volScalarField hybrid::Tmr() const
{
    volScalarField Tmr_(IOobject("TmrNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("TmrNF", dimTemperature, SMALL));
    volScalarField Tmrbf_ = baseFluidPtr_->Tmr();

    Tmr_ = Tmrbf_;

    return Tmr_;
}


const volScalarField hybrid::alpha() const
{
    volScalarField alpha_(IOobject("alphaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("alphaNF", dimless, SMALL));

    forAll(particlesComponents_, i)
    {
        alpha_ += alphasPtr_[i];
    }

    return alpha_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void hybrid::print() const
{
    volScalarField alpha_ = alpha();
    volScalarField Tm_ = Tm();
    volScalarField pv_ = pv();

    bool n1 = false;
    bool n2 = false;
    forAll(alpha_, i)
    {
        if ( T_[i] <= Tm_[i] )
        {
            n1 = true;
        }

        if ( p_[i] <= pv_[i] )
        {
            n2 = true;
        }
    }
    if ( n1 == true )
    {
        cout << boldRed << "Warning: " << reset << "Flow temperature is lower than the Fluid's melting temperature." << "\n";
    }

    if ( n2 == true )
    {
        cout << boldRed << "Warning: " << reset << "Flow pressure is lower than the Fluid's vapour pressure (i.e. saturation state)." << "\n";
    }


    boolList n3(particlesComponents_.size(), false);
    forAll(particlesComponents_, j)
    {
        volScalarField Tmnp_ = particlesProperties_[j].Tm();
        volScalarField alphanp_ = alphasPtr_[j];
        forAll(alphanp_, i)
        {
            if ( (alphanp_[i] > 0.0) & (T_[i] >= Tmnp_[i]) )
            {
                n3[j] = true;
            }
        }

        if ( n3[j] == true )
        {
            cout << boldRed << "Warning: " << reset << "Flow temperature is higher than the particle's melting temperature (i.e. particles are melting) for " << particlesProperties_[j].particleName() << "." << "\n";
        }
    }

    Info << endl;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nanoFluidModels
} // End namespace Foam

// ************************************************************************* //
