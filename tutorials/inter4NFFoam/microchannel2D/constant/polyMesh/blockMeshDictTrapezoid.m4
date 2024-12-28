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

//variables
m4_define(width, 40)
m4_define(distance, 4)
m4_define(widthu, 1)
m4_define(widthr, calc(width-distance-widthu))
m4_define(height, 1)
m4_define(heightu, 4)
m4_define(length, 1)

m4_define(numberOfCell_x, 160)
m4_define(numberOfCell_y, 16)
m4_define(numberOfCell_z, 1)
m4_define(delta_x, calc(width/numberOfCell_x))
m4_define(delta_y, calc(height/numberOfCell_y))
m4_define(delta_x, calc(length/numberOfCell_z))
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
convertToMeters 0.001;

vertices
(
  //-leftChannel
    //-Back
    (    0.0    0.0    calc(-1*length/2)    ) vlabel(Bl0)
    (    distance    0.0    calc(-1*length/2)    ) vlabel(Bl1)
    (    distance    height    calc(-1*length/2)    ) vlabel(Bl2)
    (    0.0    height    calc(-1*length/2)    ) vlabel(Bl3)

    //-Front
    (    0.0    0.0    calc(length/2)    ) vlabel(Fl0)
    (    distance    0.0    calc(length/2)    ) vlabel(Fl1)
    (    distance    height    calc(length/2)    ) vlabel(Fl2)
    (    0.0    height    calc(length/2)    ) vlabel(Fl3)

  //-rightChannel
    //-Back
    (    calc(distance+widthu)    0.0    calc(-1*length/2)    ) vlabel(Br0)
    (    width    0.0    calc(-1*length/2)    ) vlabel(Br1)
    (    width    height    calc(-1*length/2)    ) vlabel(Br2)
    (    calc(distance+widthu)    height    calc(-1*length/2)    ) vlabel(Br3)

    //-Front
    (    calc(distance+widthu)    0.0    calc(length/2)    ) vlabel(Fr0)
    (    width    0.0    calc(length/2)    ) vlabel(Fr1)
    (    width    height    calc(length/2)    ) vlabel(Fr2)
    (    calc(distance+widthu)    height    calc(length/2)    ) vlabel(Fr3)

  //-upperChannel
    //-Back
    (    calc(distance+widthu)    calc(height+heightu)    calc(-1*length/2)    ) vlabel(Bu0)
    (    distance    calc(height+heightu)    calc(-1*length/2)    ) vlabel(Bu1)

    //-Front
    (    calc(distance+widthu)    calc(height+heightu)    calc(length/2)    ) vlabel(Fu0)
    (    distance    calc(height+heightu)    calc(length/2)    ) vlabel(Fu1)
);

blocks
(
    hex (Bl0 Bl1 Bl2 Bl3 Fl0 Fl1 Fl2 Fl3) (calc(distance/delta_x) calc(height/delta_y) numberOfCell_z) simpleGrading (1 1 1)//leftChannel
    hex (Bl1 Br0 Br3 Bl2 Fl1 Fr0 Fr3 Fl2) (calc(2*widthu/delta_x) calc(height/delta_y) numberOfCell_z) simpleGrading (1 1 1)//middleChannel
    hex (Br0 Br1 Br2 Br3 Fr0 Fr1 Fr2 Fr3) (calc(widthr/delta_x) calc(height/delta_y) numberOfCell_z) simpleGrading (1 1 1)//rightChannel
    hex (Bl2 Br3 Bu0 Bu1 Fl2 Fr3 Fu0 Fu1) (calc(2*widthu/delta_x) calc(2*heightu/delta_y) numberOfCell_z) simpleGrading (1 1 1)//upperChannel
);

edges
(
);

boundary
(

    inletu
    {
        type patch;
        faces
        (
            (Bu0 Bu1 Fu1 Fu0)
        );
    }

    inlet
    {
        type patch;
        faces
        (
            (Bl0 Fl0 Fl3 Bl3)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (Br1 Br2 Fr2 Fr1)
        );
    }

    bottom
    {
        type wall;
        faces
        (
            (Bl0 Bl1 Fl1 Fl0)
            (Bl1 Br0 Fr0 Fl1)
            (Br0 Br1 Fr1 Fr0)
        );
    }

    topLeft
    {
        type wall;
        faces
        (
            (Bl2 Bl3 Fl3 Fl2)
        );
    }

    topRight
    {
        type wall;
        faces
        (
            (Br2 Br3 Fr3 Fr2)
        );
    }

    upperWalls
    {
        type wall;
        faces
        (
            (Br3 Bu0 Fu0 Fr3)
            (Bu1 Bl2 Fl2 Fu1)
        );
    }

    frontAndBack
    {
        type empty;
        faces
        (
            (Bl0 Bl3 Bl2 Bl1)
            (Bl1 Bl2 Br3 Br0)
            (Br0 Br3 Br2 Br1)
            (Fl0 Fl1 Fl2 Fl3)
            (Fl1 Fr0 Fr3 Fl2)
            (Fr0 Fr1 Fr2 Fr3)
        );
    }

    frontAndBacku
    {
        type empty;
        faces
        (
            (Bl2 Bu1 Bu0 Br3)
            (Fu0 Fu1 Fl2 Fr3)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
