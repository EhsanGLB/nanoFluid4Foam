
if(typeFluid[i] == "forcedConvection")
{
    #include "forcedConvection/TEqnFC.H"

    if (nanoFluidModelPtr.nanoFluidModelType())
    {
        #include "forcedConvection/alphaEqnFC.H"
    }
}


if(typeFluid[i] == "naturalConvection")
{
    if(PCMProperties[i].headerOk())
    {
        #include "PCM/TEqnPCM.H"

        if (nanoFluidModelPtr.nanoFluidModelType())
        {
            #include "PCM/alphaEqnPCM.H"
        }
    }

    else
    {
        #include "naturalConvection/TEqnNC.H"

        if (nanoFluidModelPtr.nanoFluidModelType())
        {
            #include "naturalConvection/alphaEqnNC.H"
        }
    }
}


