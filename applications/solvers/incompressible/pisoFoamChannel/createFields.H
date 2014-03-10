    // Create and read in initial velocity field.
    Info << "Creating/reading velocity field, U..." << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Compute the fluxes.
    Info << "Creating/calculating face flux field, phi..." << endl; 
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );

    // Create and read in initial pressure field.
    Info << "Creating/reading pressure field, p..." << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Create access to transport properties dictionary.
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    // Create and read in uniform viscosity.
    Info << "Reading viscosity..." << endl;
    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    // Create and read in the bulk channel velocity.
    Info << "Reading bulk channel velocity..." << endl;
    dimensionedVector Ubar
    (
        transportProperties.lookup("Ubar")
    );
