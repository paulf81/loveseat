    dimensionedVector gradP
    (
        "gradP",
        dimensionSet(0, 1, -2, 0, 0),
        vector::zero
    );

    IFstream gradPFile
    (
        runTime.path()/runTime.timeName()/"uniform"/"gradP.raw"
    );

    if(gradPFile.good())
    {
        gradPFile >> gradP;
    //  gradP = (gradP & flowDirectionVector) * flowDirectionVector;
        Info << "Reading average pressure gradient..." << endl;
    }
    else
    {
        Info<< "Initializing with 0 pressure gradient..." << endl;
    };
