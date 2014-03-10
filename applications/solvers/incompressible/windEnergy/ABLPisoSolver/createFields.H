    // Create and read in initial velocity field
    Info<< "Creating/Reading velocity field, U..." << endl;
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


    // Create initial average velocity field
    Info<< "Creating mean velocity field, Umean..." << endl;
    volVectorField Umean
    (
        IOobject
        (
            "Umean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        U
    );


    // Create initial fluctuating velocity field
    Info<< "Creating fluctuating velocity field, Uprime..." << endl;
    volVectorField Uprime
    (
        IOobject
        (
            "Uprime",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U
    );


    // Create divergence of velocity field
    Info<< "Creating divergence of velocity field..." << endl;
    volScalarField divU
    (
        IOobject
        (
            "divU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("divU",dimensionSet(0, 0, -1, 0, 0, 0, 0),0.0)
    );


    // Create and read in initial temperature field
    Info<< "Reading temperature field, T..." << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Create initial average temperature field
    Info<< "Creating mean temperature field, Tmean..." << endl;
    volScalarField Tmean
    (
        IOobject
        (
            "Tmean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        T
    );

    // Create initial fluctuating temperature field
    Info<< "Creating fluctuating temperature field, Tprime..." << endl;
    volScalarField Tprime
    (
        IOobject
        (
            "Tprime",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        T
    );


    // Create and read in initial pressure field
    Info<< "Creating/Reading deviation from hydrostatic pressure field, pd..." << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // Pressure reference information
    Info<< "Reading pressure field reference information..." << endl;
    bool pdRefOn = mesh.solutionDict().subDict("options").lookupOrDefault<bool>("pdRefOn", true);
    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell
    (
        pd,
        mesh.solutionDict().subDict("options"),
        pdRefCell,
        pdRefValue
    );

    // Compute the fluxes
    Info<< "Creating/Calculating face flux field, phi..." << endl;
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

    // Create divergence of velocity flux
    Info<< "Creating divergence of velocity flux..." << endl;
    volScalarField divPhi
    (
        IOobject
        (
            "divPhi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("divPhi",dimensionSet(0, 0, -1, 0, 0, 0, 0),0.0)
    );


    // Get transport properties
    Info << "Reading transport properties..." << endl;
    #include "readTransportProperties.H"

    
    // Get atmospheric boundary layer properties
    Info << "Reading atmospheric boundary layer properties..." << endl;
    #include "readABLProperties.H"

    // Create kinematic density field for buoyancy force
    Info<< "Creating the kinematic density field, rhok..." << endl;
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - (T - TRef)/TRef
    );


    // Create Coriolis force vector
    Info<< "Creating the Coriolis force vector, fCoriolis..." << endl;
    volVectorField fCoriolis
    (
        IOobject
        (
            "fCoriolis",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
	mesh,
	dimensionedVector("fCoriolis",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
    );


    // Create the SGS energy
    Info<< "Creating the turbulent kinetic energy field, kLES..." << endl;
    surfaceScalarField kLES
    (
        IOobject
        (
            "kLES",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("kLES",dimensionSet(0, 2, -2, 0, 0, 0, 0),0.0)
    );


    // Create the momentum diffusivity field
    Info<< "Creating the turbulent viscosity field, nuLES..." << endl;
    surfaceScalarField nuLES
    (
        IOobject
        (
            "nuLES",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("nuLES",dimensionSet(0, 2, -1, 0, 0, 0, 0),0.0)
    );


    // Create turbulent Prandtl number field
    Info<< "Creating the turbulent Prandtl number field, Prt..." << endl;
    surfaceScalarField Prt
    (
        IOobject
        (
            "Prt",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Prt",dimensionSet(0,0,0,0,0,0,0),1.0/3.0)
    );


    // Create thermal diffusivity field
    Info<< "Creating the turbulent thermal diffusivity field, kappaLES..." << endl;
    surfaceScalarField kappaLES
    (
        IOobject
        (
            "kappaLES",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	nuLES/Prt
    );


    // Create the deviatoric part of the SGS momentum flux field
    Info<< "Creating the deviatoric SGS stresses field, devR..." << endl;
    surfaceSymmTensorField devR
    (
        IOobject
        (
            "devR",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedSymmTensor("devR",dimensionSet(0, 2, -2, 0, 0, 0, 0),symmTensor::zero)
    );


    // Create the divergence of the deviatoric part of the SGS momentum flux field
    Info<< "Creating the divergence of the deviatoric SGS stresses field, divDevR..." << endl;
    volVectorField divDevR
    (
        IOobject
        (
            "divDevR",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("divDevR",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
    );


    // Create the SGS temperature flux field
    Info<< "Creating the SGS temperature flux field, q..." << endl;
    surfaceVectorField q
    (
        IOobject
        (
            "q",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("q",dimensionSet(0, 1, -1, 1, 0, 0, 0),vector::zero)
    );


    // Create the divergence of the SGS temperature flux field
    Info<< "Creating the divergence of the SGS temperature flux field, divq..." << endl;
    volScalarField divq
    (
        IOobject
        (
            "divq",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("divq",dimensionSet(0, 0, -1, 1, 0, 0, 0),0.0)
    );
