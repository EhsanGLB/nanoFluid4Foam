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

//**************************************************************************# Ceramic material #**************************************************************************//
AL2O3    newMaterialConst    coeffs    101.96    1.0e-23    1.0e-23    1.0e-23    2327.0    1.0    1160000.0    3970.0    40.0    765.0    0.85e-5    5.0e23    45.0e-3    1.0e5;


//**************************************************************************# Organic and Carbon-based material #***********************************************************//


//**************************************************************************# Polymeric material #**************************************************************************//

//**************************************************************************# Inorganic material #**************************************************************************//
water    newMaterialConst    coeffs    18.015    647.13    2.2055e+7    0.05595    273.15    1.0    333606.0    993    0.628    4178.0    36.2e-5    2.0e+7    696e-6    1.0e6;//25 C


