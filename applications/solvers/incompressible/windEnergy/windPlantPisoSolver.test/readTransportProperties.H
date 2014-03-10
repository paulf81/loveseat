    // Molecular viscosity [m^2/s]
    dimensionedScalar nu(laminarTransport.lookup("nu"));

    // Reference temperature [K]
    dimensionedScalar TRef(laminarTransport.lookup("TRef"));

    // Von Karman constant for use in Monin-Obuhkov similarity
    scalar kappa(laminarTransport.lookupOrDefault<scalar>("kappa", 0.41));

    // Constants for Monin-Obhukov universal functions
    scalar betaM(laminarTransport.lookupOrDefault<scalar>("betaM", 16.0));
    scalar gammaM(laminarTransport.lookupOrDefault<scalar>("gammaM", 5.0));

    // Surface roughness height [m]
    dimensionedScalar z0(laminarTransport.lookup("z0"));

    // Surface mean temperature flux [K-m/s]
    dimensionedScalar q0(laminarTransport.lookup("q0"));

    // Surface stress model (options are "Schummann", "Moeng", or "MoengCorrected")
    word surfaceStressModel(laminarTransport.lookupOrDefault<word>("surfaceStressModel","Moeng"));

    // Wei-Brasseur surface stress model correction beta parameter
    scalar betaSurfaceStress(laminarTransport.lookupOrDefault<scalar>("betaSurfaceStress", 0.15));
