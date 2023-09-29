rho = nanoFluidModelPtr->rho();
kappa = nanoFluidModelPtr->kappa();
Cp = nanoFluidModelPtr->Cp();
mu = nanoFluidModelPtr->mu();
beta = nanoFluidModelPtr->beta();
nanoFluidModelPtr->print();

forAll(moltensPtr, i)
{
    moltensPtr[i] = particlesProperties[i].molten();
}

