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

#include "NFDimensionSet.H"//-nanoFluid4Foam
#include "nanoFluidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nanoFluidModel, 0);
    defineRunTimeSelectionTable(nanoFluidModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nanoFluidModel::nanoFluidModel
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T,
    const PtrList<volScalarField>& alphasPtr
)
:
    nanoFluidPropertiesDict_(nanoFluidPropertiesDict),
    U_(U),
    p_(p),
    T_(T),
    alphasPtr_(alphasPtr),
    nanoFluidModelList_(nanoFluidPropertiesDict.lookup("nanoFluidModel")),
    nanoFluidModelType_(nanoFluidModelList_[1]),
    baseFluidPtr_(baseFluid::New(nanoFluidPropertiesDict, U, p, T)),
    particlesComponents_(nanoFluidPropertiesDict.lookup("particle")),
    particlesProperties_(particlesComponents_.size()),
    DBs_(particlesComponents_.size()),
    DTs_(particlesComponents_.size())
{
    forAll(particlesProperties_, i)
    {
            particlesProperties_.set(i, particleModel::New(particlesComponents_[i], nanoFluidPropertiesDict, U, p, T).ptr() );
    }


    autoPtr<brownianDiffusivityModel> brownianDiffusivityModelPtr_(brownianDiffusivityModel::New(nanoFluidPropertiesDict_, U_, p_, T_, alphasPtr_, baseFluidPtr_, particlesProperties_).ptr() );
    brownianDiffusivityModel& brownianDiffusivityModelObj = brownianDiffusivityModelPtr_();
    const PtrList<volScalarField>& DBss_ = brownianDiffusivityModelObj.DBs();
    forAll(DBs_, i)
    {
            DBs_.set(i, DBss_[i] );
    }


    autoPtr<thermophoresisDiffusivityModel> thermophoresisDiffusivityModelPtr_(thermophoresisDiffusivityModel::New(nanoFluidPropertiesDict_, U_, p_, T_, alphasPtr_, baseFluidPtr_, particlesProperties_).ptr() );
    thermophoresisDiffusivityModel& thermophoresisDiffusivityModelObj = thermophoresisDiffusivityModelPtr_();
    const PtrList<volScalarField>& DTss_ = thermophoresisDiffusivityModelObj.DTs();
    forAll(DTs_, i)
    {
            DTs_.set(i, DTss_[i] );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*
inline const PtrList<volScalarField>& nanoFluidModel::DBs() const
{
    return DBs_;
}
*/

bool nanoFluidModel::read(const dictionary& nanoFluidPropertiesDict)
{
    nanoFluidPropertiesDict_ = nanoFluidPropertiesDict;

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
