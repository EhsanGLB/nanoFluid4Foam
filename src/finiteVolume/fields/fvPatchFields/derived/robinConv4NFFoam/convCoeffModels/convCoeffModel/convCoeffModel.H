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

Namespace
    Foam::convCoeffModels

Description
    A namespace for various incompressible convCoeffModel implementations.

Class
    Foam::convCoeffModel

Description
    An abstract base class for incompressible convCoeffModels

SourceFiles
    convCoeffModel.C
    newconvCoeffModel.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef convCoeffModel_H
#define convCoeffModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class convCoeffModel Declaration
\*---------------------------------------------------------------------------*/

class convCoeffModel
{
protected:

    // Protected data
        dictionary dict_;
        word convCoeffModelType_;
        scalar To_;
        scalarList flowProp_;
        scalar rhoc_;
        scalar kappac_;
        scalar Cpc_;
        scalar muc_;
        scalar betac_;



    // Private Member Functions

        //- Disallow copy construct
        convCoeffModel(const convCoeffModel&);

        //- Disallow default bitwise assignment
        void operator=(const convCoeffModel&);


public:

    //- Runtime type information
    TypeName("convCoeffModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            convCoeffModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Return a reference to the selected baseFluid model
        static autoPtr<convCoeffModel> New
        (
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        convCoeffModel
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~convCoeffModel()
        {}


    // Member Functions

        //- Return the dictBC
        const dictionary& dictBC() const
        {
            return dict_;
        }

        //- density [kg/m^3]
        scalar rho() const;

        //- thermal conductivity [W/(m K)]
        scalar kappa() const;

        //- heat capacity [J/(kg K)]
        scalar Cp() const;

        //- dynamic viscosity [Pa s]
        scalar mu() const;

        //- thermal expansion [1/K]
        scalar beta() const;

        //- kinematic viscosity [m^2/s]
        scalar nu() const;

        //- thermal diffusivity [m^2/s]
        scalar alpha() const;

        //- Prandtl number [-]
        scalar Pr() const;



        //- Reynolds number [-]
        //virtual scalar Re() const = 0;

        //- Rayleigh number [-]
        //virtual scalar Ra() const = 0;

        //- Rayleigh number [-]
        //virtual scalar Gr() const = 0;

        //- Return the average convective heat transfer coefficient [W/m^2 K]
        virtual scalar hoAvg(scalar avgT) const = 0;

        //- Write
        virtual void write(Ostream& os) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
