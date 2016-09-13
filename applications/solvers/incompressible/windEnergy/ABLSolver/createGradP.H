    // Create the driving pressure gradient field.    
    dimensionedVector gradP
    (
        "gradP",
        dimensionSet(0, 1, -2, 0, 0, 0, 0),
        vector::zero
    );

    // Create the driving pressure gradient time history file.
    if (Pstream::master() && !isDir(postProcessingDir/"gradPHistory"))
    {
        mkDir(postProcessingDir/"gradPHistory");
    }
    if (Pstream::master() && !isDir(postProcessingDir/"gradPHistory"/runTime.timeName()))
    {
        mkDir(postProcessingDir/"gradPHistory"/runTime.timeName());
    }

    autoPtr<OFstream> gradPHistoryFile(NULL);
    if (Pstream::master())
    {
        gradPHistoryFile.reset(new OFstream(postProcessingDir/"gradPHistory"/runTime.timeName()/"gradPHistory"));
        gradPHistoryFile() << "Time (s)" << " " << "pressure gradient (m/s^2)" << " "
                                                << "Desired Wind Speed (m/s)" << " " 
                                                << "Desired Wind Direction (deg)" << " "
                                                << "Desired Wind Vector (m/s)" << " "
                                                << "Uncorrected Wind Vector (m/s)" << endl;
    }
