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
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("UMean",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
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
        mesh,
	dimensionedVector("Uprime",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    );


    // Create the resolved Reynolds stress tensor field
    Info<< "Creating the resolved Reynolds stress field, UprimeUprimeMean..." << endl;
    volSymmTensorField UprimeUprimeMean
    (
        IOobject
        (
            "UprimeUprimeMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedSymmTensor("UprimeUprimeMean",dimensionSet(0, 2, -2, 0, 0, 0, 0),symmTensor::zero)
    );


    // Create the resolved turbulent kinetic energy field
    Info<< "Creating the resolved turbulent kinetic energy field, k..." << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0.5 * tr(UprimeUprimeMean)
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
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("Tmean",dimensionSet(0, 0, 0, 1, 0, 0, 0),0.0)
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
        mesh,
	dimensionedScalar("Tprime",dimensionSet(0, 0, 0, 1, 0, 0, 0),0.0)
    );


    // Create the temperature-velocity correlation field
    Info<< "Creating the temperature-velocity correlation field, TprimeUprimeMean..." << endl;
    volVectorField TprimeUprimeMean
    (
        IOobject
        (
            "TprimeUprimeMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("TprimeUprimeMean",dimensionSet(0, 1, -1, 1, 0, 0, 0),vector::zero)
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


    // Get transport properties
    Info << "Reading transport properties..." << endl;
    #include "readTransportProperties.H"

    
    // Get atmospheric boundary layer properties
    Info << "Reading atmospheric boundary layer properties..." << endl;
    #include "readABLProperties.H"

    // Create an object of the horizontalWindTurbineArray class if there
    // is to be a turbine array
    turbineModels::horizontalAxisWindTurbinesALMold turbines(U);
    
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


    //   Calculate the LES filter width (assumes the mesh is uniform spacing throughout).
    Info<< "Calculating LES filter widths, deltaLES..." << endl;
    #include "setDeltaLES.H"


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
            IOobject::NO_WRITE
        ),
	nuLES/(1.0/3.0)
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
