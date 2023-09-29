

    fvScalarMatrix alphamEqn
    ( 
          fvm::div(phi, alpham)
        - fvm::laplacian(Dbm, alpham)
        + fvm::laplacian(((W*b*alpham)/(alphasPtr[0] - CRef)), alphasPtr[0])
    );

    alphamEqn.relax();
    alphamEqn.solve();
//******************mp Equation ***********************************************************
