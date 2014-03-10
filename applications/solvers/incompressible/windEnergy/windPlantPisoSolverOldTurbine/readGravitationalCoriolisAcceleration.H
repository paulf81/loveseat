    // Read gravitational acceleration vector
    Info << "Reading gravitational acceleration..." << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    // Read rotation rate vector of the planetary body
    Info << "Reading planetary body rotation rate..." << endl;
    uniformDimensionedVectorField Omega
    (
        IOobject
        (
            "Omega",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
