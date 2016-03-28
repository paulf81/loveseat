    // Create the source term, which lumps together the driving pressure gradient and mesoscale source terms
    // for advection of momentum and potential temperature.

    // Create the momentum sourc term.
    Info << "Creating the momentum source term, SourceU..." << endl;
    volVectorField SourceU
    (
        IOobject
        (
            "SourceU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("SourceU",dimVelocity/dimTime,vector::zero)
    );

    Info << "Creating the potential temperature source term, SourceT..." << endl;
    volScalarField SourceT
    (
        IOobject
        (
            "SourceT",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SourceT",dimTemperature/dimTime,0.0)
    );



    // Create the momentum source term time history files.
    if (Pstream::master() && !isDir(postProcessingDir/"SourceUHistory"))
    {
        mkDir(postProcessingDir/"SourceUHistory");
    }
    if (Pstream::master() && !isDir(postProcessingDir/"SourceUHistory"/runTime.timeName()))
    {
        mkDir(postProcessingDir/"SourceUHistory"/runTime.timeName());
    }

    autoPtr<OFstream> sourceUHistoryFile(NULL);
    if (Pstream::master())
    {
        sourceUHistoryFile.reset(new OFstream(postProcessingDir/"sourceUHistory"/runTime.timeName()/"sourceUHistory"));
        sourceUHistoryFile() << "Time (s)" << " " << "source term (m/s^2)" << endl;
    }


    // Create the potential temperature source term time history files.
    if (Pstream::master() && !isDir(postProcessingDir/"SourceTHistory"))
    {
        mkDir(postProcessingDir/"SourceTHistory");
    }
    if (Pstream::master() && !isDir(postProcessingDir/"SourceTHistory"/runTime.timeName()))
    {
        mkDir(postProcessingDir/"SourceTHistory"/runTime.timeName());
    }

    autoPtr<OFstream> sourceTHistoryFile(NULL);
    if (Pstream::master())
    {
        sourceTHistoryFile.reset(new OFstream(postProcessingDir/"sourceTHistory"/runTime.timeName()/"sourceTHistory"));
        sourceTHistoryFile() << "Time (s)" << " " << "source term (K/s)" << endl;
    }
