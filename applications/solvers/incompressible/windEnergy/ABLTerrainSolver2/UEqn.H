    // Solve the momentum equation

    #include "computeCoriolisForce.H"

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)                             // time derivative
      + fvm::div(phi, U)                        // convection
      + turbulence->divDevReff(U)               // momentum flux
      + fvc::div(Rwall)
      - fCoriolis                               // Coriolis force
      + gradP                                   // driving pressure gradient
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
            ==
            fvc::reconstruct
            (
                - fvc::snGrad(p_rgh)*mesh.magSf()    // pressure gradient
                +(fvc::interpolate(rhok - 1.0)*(g & mesh.Sf()))  // Boussinesq buoyancy term
            )
        );
    }
