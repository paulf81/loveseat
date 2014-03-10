    // Create the vorticity field.
    Info<< "Creating the vorticity Field, omega..." << endl;
    volVectorField omega
    (
        IOobject
        (
            "omega",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(U)
    );


    // Create the second-invariant field.
    Info<< "Creating the second-invariant of the velocity-gradient tensor field, Q..." << endl;
    volScalarField Q
    (
        IOobject
        (
            "Q",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*(sqr(tr(fvc::grad(U))) - tr(((fvc::grad(U))&(fvc::grad(U)))))
    );

