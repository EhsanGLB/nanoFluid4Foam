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

Class
    Foam::LNH4NFFoam

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates the LNH4NFFoam, the curl of the velocity.

SourceFiles
    LNH4NFFoam.C
    IOLNH4NFFoam.H

\*---------------------------------------------------------------------------*/

#ifndef LNH4NFFoam_H
#define LNH4NFFoam_H

#include "volFieldsFwd.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class LNH4NFFoam Declaration
\*---------------------------------------------------------------------------*/

class LNH4NFFoam
{
    // Private data

        //- Name of this set of LNH4NFFoam objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of velocity field, default is "U"
        word UName_;

        //- Name of LNH4NFFoam field
        word outputName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        LNH4NFFoam(const LNH4NFFoam&);

        //- Disallow default bitwise assignment
        void operator=(const LNH4NFFoam&);


public:

    //- Runtime type information
    TypeName("LNH4NFFoam");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        LNH4NFFoam
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~LNH4NFFoam();


    // Member Functions

        //- Return name of the set of LNH4NFFoam
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the LNH4NFFoam data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the LNH4NFFoam and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
