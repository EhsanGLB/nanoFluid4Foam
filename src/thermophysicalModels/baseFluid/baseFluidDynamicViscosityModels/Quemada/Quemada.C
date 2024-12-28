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
#include "Quemada.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace baseFluidDynamicViscosityModels
{
    defineTypeNameAndDebug(Quemada, 0);
    addToRunTimeSelectionTable(baseFluidDynamicViscosityModel, Quemada, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Quemada::Quemada
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    baseFluidDynamicViscosityModel(nanoFluidPropertiesDict, U, p, T),
    QuemadaCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    a0_(QuemadaCoeffs_.lookup("a0")),
    a1_(QuemadaCoeffs_.lookup("a1")),
    b0_(QuemadaCoeffs_.lookup("b0")),
    b1_(QuemadaCoeffs_.lookup("b1")),
    b2_(QuemadaCoeffs_.lookup("b2")),
    b3_(QuemadaCoeffs_.lookup("b3")),
    c0_(QuemadaCoeffs_.lookup("c0")),
    c1_(QuemadaCoeffs_.lookup("c1")),
    c2_(QuemadaCoeffs_.lookup("c2")),
    c3_(QuemadaCoeffs_.lookup("c3")),
    gammaC0_(QuemadaCoeffs_.lookup("gammaC0")),
    muPlasma_(QuemadaCoeffs_.lookup("muPlasma"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Quemada::mu() const
{
    volScalarField mu_(IOobject("muBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muBF", dimDynamicViscosity, SMALL));

    const volScalarField& H= U_.mesh().lookupObject<volScalarField>("H");

    mu_ = muPlasma_ *
        pow
        (
         (
          1.0-0.5*
          // K: Quemada Parameter
          (
           (
            (a0_+2/(a1_+H))                                   // k0
            +
            exp(b0_+b1_*H+b2_*pow(H,2.0)+b3_*pow(H,3.0))        // kInf
            *sqrt
            (
             strainRate()/
             (
              exp(c0_+c1_*H+c2_*pow(H,2.0)+c3_*pow(H,3.0))
              *gammaC0_                                        // gammaC
             )
            )
           )
           /
           (
            1+sqrt
            (
             strainRate()/
             (
              exp(c0_+c1_*H+c2_*pow(H,2.0)+c3_*pow(H,3.0))
              *gammaC0_                                        // gammaC
             )
            )
           )
           )
           // end K
           *H
           ),-2
           )
           ;

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace baseFluidDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
