//------------------------------- nanoFluid4Foam project -------------------------------//
//Author
    //Ehsan Golab, SUT. All rights reserved.
    //Ehsan1996Golab@gmail.com

//--------------------------------------------------------------------------------------//

//- Templates
//- A material with const properties
//{
//Name    newMaterialConst    coeffs    molecularWeight[kg/kmol]    criticalTemperature[K]    criticalPressure[Pa]    criticalVolume[m^3/mol]    meltingPoint[K]    meltingRange[K]    meltingEnthalpy[J/kg]    density[kg/m^3]    ThermalConductivity[W/(m.K)]    specificHeatCapacity[J/(kg.K)]    thermalExpansion[1/K]    electricalResistivity[ohm.m]    dynamicViscosity[Pa.s]    vapourPressure[Pa];
//}

//- A material with non-const properties
//NameDict
//{
//    newMaterialFuncT    newMaterialFuncT;
//    constCoeffs    molecularWeight[kg/kmol]    criticalTemperature[K]    criticalPressure[Pa]    criticalVolume[m^3/mol]    meltingPoint[K]    meltingRange[K]    meltingEnthalpy[J/kg]
//    rhoCoeffs      a0    a1    a2    a3    a4    a5    a6    b    c;
//    KCoeffs        a0    a1    a2    a3    a4    a5    a6    b    c;
//    cpCoeffs       a0    a1    a2    a3    a4    a5    a6    b    c;
//    betaCoeffs     a0    a1    a2    a3    a4    a5    a6    b    c;
//    rhoRCoeffs     a0    a1    a2    a3    a4    a5    a6    b    c;
//    muCoeffs       a0    a1    a2    a3    a4    a5    a6    b    c;
//    pvCoeffs       a0    a1    a2    a3    a4    a5    a6    b    c;
//}
//- f(T) = ( a0 + a1*T^1 + a2*T^2 + a3*T^3 + a4*T^4 + a5*T^5 + a6*T^6 ) * exp( b * T ) + c;


//**************************************************************************# Metallic material #**************************************************************************//
Cu    newMaterialConst    coeffs    63.54    1.0e-23    1.0e-23    1.0e-23    1356.0    2.0    205382.0    8900.0    398.0    386.0    17.6e-6    17.1e-9    4.25e-3    1.0e5;
CuFuncTDict
{
    newMaterialFuncT    newMaterialFuncT;
    constCoeffs    63.54    1.0e-23    1.0e-23    1.0e-23    1356.0    2.0    205382.0;
    rhoCoeffs      8900   0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    KCoeffs        441.6    -0.17119    -1.5446e-04    -7.2917e-8    0.0    0.0    0.0    0.0    0.0;
    cpCoeffs       285.8    0.44631    -5.2054e-04    2.3958e-07    0.0    0.0    0.0    0.0    0.0;
    betaCoeffs     1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    rhoRCoeffs     1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    muCoeffs       1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    pvCoeffs       1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
}


//**************************************************************************# Ceramic material #**************************************************************************//
AL2O3    newMaterialConst    coeffs    101.96    1.0e-23    1.0e-23    1.0e-23    2327.0    1.0    1160000.0    3970.0    25.0    703.0    7.7e-6    5.0e23    45.0e-3    1.0e5;
AL2O3FuncTDict
{
    newMaterialFuncT    newMaterialFuncT;
    constCoeffs    101.96    1.0e-23    1.0e-23    1.0e-23    2327.0    1.0    1160000.0;
    rhoCoeffs      3970   0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    KCoeffs        148.14    -0.56883    9.794e-04    -8.0417e-7    2.5595e-10    0.0    0.0    0.0    0.0;
    cpCoeffs       531.43    7.135    -0.011923    9.3125e-06    -2.7679e-09    0.0    0.0    0.0    0.0;
    betaCoeffs     1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    rhoRCoeffs     1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    muCoeffs       1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    pvCoeffs       1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
}


//**************************************************************************# Organic and Carbon-based material #***********************************************************//
MWCNT    newMaterialConst    coeffs    12.01    1.0e-23    1.0e-23    1.0e-23    3823.0    2.0    1e-10    1600.0    3000.0    796.0    8.5e-6    1.0e-7    4.25e-3    1.0e5;


//**************************************************************************# Polymeric material #**************************************************************************//


//**************************************************************************# Inorganic material #**************************************************************************//
therminolVP1    newMaterialConst    coeffs    166.0    722.0    33.1e5    0.003058104    285.0    5.0    4700    1056.0    0.1353    1575.0    1.0e-23    1.0e-23    3.28e-3    4.0;//300 K
//therminolVP1    newMaterialConst    coeffs    166.0    722.0    33.1e5    0.003058104    285.0    5.0    4700    974.0    0.1239    1858.0    1.0e-23    1.0e-23    0.707e-3    2000.0;//400 K
//therminolVP1    newMaterialConst    coeffs    166.0    722.0    33.1e5    0.003058104    285.0    5.0    4700    886.0    0.1089    2128.0    1.0e-23    1.0e-23    0.324e-3    53600;//500 K

therminolVP1FuncTDict
{
//285 < T < 693 K in liquid phase
    newMaterialFuncT    newMaterialFuncT;
    constCoeffs    166.0    722.0    33.1e5    0.003058104    285.0    5.0    4700;
    rhoCoeffs      1.4386e+03    -1.8711e+00    2.737e-03    -2.3793e-06    0.0    0.0    0.0    0.0    0.0;
    KCoeffs        0.14644    2.0353e-05    -1.9367e-07    1.0614e-11    0.0    0.0    0.0    0.0    0.0;
    cpCoeffs       2.125e+03    -11.017    0.049862    -7.7663e-05    4.394e-08    0.0    0.0    0.0    0.0;
    betaCoeffs     1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    rhoRCoeffs     1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    muCoeffs       4.15    0.0    0.0    0.0    0.0    0.0    0.0    -2.363e-02    2.5847e-04;
    pvCoeffs       1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
}
