/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     4.1                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//------------------------------- nanoFluid4Foam project -------------------------------//
//Author
    //Ehsan Golab, SUT. All rights reserved.
    //Ehsan1996Golab@gmail.com

//--------------------------------------------------------------------------------------//

//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])

m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

m4_define(pi, 3.14159265358979323844)
m4_define(rad, [calc($1*pi/180.0)])
m4_define(startAngleOffset, rad(45.0))
m4_define(rOut, calc(tan(60.0*pi/180.0)))



//dimensions of internal rectangular [x, y, z]->[W, H, L]
//width and height of rectangular
m4_define(intRecW, 0.3)
m4_define(intRecH, 0.25)
m4_define(intRecL, 0.1)
m4_define(intRecTheta, calc(atan(2*intRecH/intRecW)*180/pi))

//Number of cell
m4_define(intRecWN, 20)
m4_define(intRecHN, 20)
m4_define(intRecLN, 1)

//dimensions of internal cylinder
m4_define(intCyR, 0.5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
convertToMeters 1.0;

vertices
(
  //-rectangular
    //-Back of rectangular
    (    calc(intRecW/2)    0    calc(-1*intRecL/2)    ) vlabel(intRecB0)
    (    calc(intRecW/2)    calc(-1*intRecH)    calc(-1*intRecL/2)    ) vlabel(intRecB1)
    (    calc(-1*intRecW/2)    calc(-1*intRecH)    calc(-1*intRecL/2)    ) vlabel(intRecB2)
    (    calc(-1*intRecW/2)    0    calc(-1*intRecL/2)    ) vlabel(intRecB3)
    //-Front of rectangular
    (    calc(intRecW/2)    0    calc(intRecL/2)    ) vlabel(intRecF0)
    (    calc(intRecW/2)    calc(-1*intRecH)    calc(intRecL/2)    ) vlabel(intRecF1)
    (    calc(-1*intRecW/2)    calc(-1*intRecH)    calc(intRecL/2)    ) vlabel(intRecF2)
    (    calc(-1*intRecW/2)    0    calc(intRecL/2)    ) vlabel(intRecF3)

  //-cylinder
    //-Back of cylinder
    (    calc(intCyR)    0    calc(-1*intRecL/2)    ) vlabel(intCyB0)
    (    calc(0.5*cos(rad(intRecTheta)))    calc(-1*0.5*sin(rad(intRecTheta)))    calc(-1*intRecL/2)    ) vlabel(intCyB1)
    (    calc(-1*0.5*cos(rad(intRecTheta)))    calc(-1*0.5*sin(rad(intRecTheta)))    calc(-1*intRecL/2)    ) vlabel(intCyB2)
    (    calc(-1*intCyR)    0    calc(-1*intRecL/2)    ) vlabel(intCyB3)
    //-Front of cylinder
    (    calc(intCyR)    0    calc(intRecL/2)    ) vlabel(intCyF0)
    (    calc(0.5*cos(rad(intRecTheta)))    calc(-1*0.5*sin(rad(intRecTheta)))    calc(intRecL/2)    ) vlabel(intCyF1)
    (    calc(-1*0.5*cos(rad(intRecTheta)))    calc(-1*0.5*sin(rad(intRecTheta)))    calc(intRecL/2)    ) vlabel(intCyF2)
    (    calc(-1*intCyR)    0    calc(intRecL/2)    ) vlabel(intCyF3)
);

blocks
(
    hex (intRecB0 intRecB1 intRecB2 intRecB3 intRecF0 intRecF1 intRecF2 intRecF3) (intRecWN intRecHN intRecLN) simpleGrading (1 1 1)
    hex (intRecB0 intCyB0 intCyB1 intRecB1 intRecF0 intCyF0 intCyF1 intRecF1) (intRecWN intRecHN intRecLN) simpleGrading (1 1 1)
    hex (intRecB1 intCyB1 intCyB2 intRecB2 intRecF1 intCyF1 intCyF2 intRecF2) (intRecWN intRecHN intRecLN) simpleGrading (1 1 1)
    hex (intRecB2 intCyB2 intCyB3 intRecB3 intRecF2 intCyF2 intCyF3 intRecF3) (intRecWN intRecHN intRecLN) simpleGrading (1 1 1)
);

edges
(
  //-Internal cylinder
    //-Back of internal cylinder
    arc intCyB0 intCyB1 (    calc(0.5*cos(rad(intRecTheta/2)))    calc(-1*0.5*sin(rad(intRecTheta/2)))    calc(-1*intRecL/2)    )
    arc intCyB1 intCyB2 (    0    calc(-1*intCyR)    calc(-1*intRecL/2)    )
    arc intCyB2 intCyB3 (    calc(-1*0.5*cos(rad(intRecTheta/2)))    calc(-1*0.5*sin(rad(intRecTheta/2)))    calc(-1*intRecL/2)    )

    //-Front of internal cylinder
    arc intCyF0 intCyF1 (    calc(0.5*cos(rad(intRecTheta/2)))    calc(-1*0.5*sin(rad(intRecTheta/2)))    calc(intRecL/2)    )
    arc intCyF1 intCyF2 (    0    calc(-1*intCyR)    calc(intRecL/2)    )
    arc intCyF2 intCyF3 (    calc(-1*0.5*cos(rad(intRecTheta/2)))    calc(-1*0.5*sin(rad(intRecTheta/2)))    calc(intRecL/2)    )
);

boundary
(
    horizontalWall
    {
        type wall;
        faces
        (
          //-Internal rectangular
            (intRecB0 intRecB3 intRecF3 intRecF0)
          //-Internal cylinder
            (intCyB0 intRecB0 intRecF0 intCyF0)
            (intCyB3 intRecB3 intRecF3 intCyF3)
        );
    }

    semiCylindericalWall
    {
        type wall;
        faces
        (
          //-Internal cylinder
            (intCyB1 intCyB0 intCyF0 intCyF1)
            (intCyB2 intCyB1 intCyF1 intCyF2)
            (intCyB3 intCyB2 intCyF2 intCyF3)
        );
    }

    back
    {
        type empty;
        faces
        (
          //-Internal rectangular
            (intRecB0 intRecB1 intRecB2 intRecB3)
          //-Internal cylinder
            (intRecB0 intCyB0 intCyB1 intRecB1)
            (intRecB1 intCyB1 intCyB2 intRecB2)
            (intRecB2 intCyB2 intCyB3 intRecB3)
        );
    }

    front
    {
        type empty;
        faces
        (
          //-Internal rectangular
            (intRecF0 intRecF1 intRecF2 intRecF3)
          //-Internal cylinder
            (intRecF0 intCyF0 intCyF1 intRecF1)
            (intRecF1 intCyF1 intCyF2 intRecF2)
            (intRecF2 intCyF2 intCyF3 intRecF3)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
