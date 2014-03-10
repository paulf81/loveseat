    // Pressure reference information.
    bool pRefOn = mesh.solutionDict().subDict("options").lookupOrDefault<bool>("pRefOn", true);
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        mesh.solutionDict().subDict("options"),
        pRefCell,
        pRefValue
    );
    
    // Read information about statistics gathering frequency and the direction of homogeneity.
    bool horizontalStatisticsOn = mesh.solutionDict().subDict("options").lookupOrDefault<bool>("horizontalStatisticsOn", false);
    int horizontalStatisticsFreq = mesh.solutionDict().subDict("options").lookupOrDefault<int>("horizontalStatisticsFrequency", 0);
    scalar horizontalStatisticsStartTime = mesh.solutionDict().subDict("options").lookupOrDefault<scalar>("horizontalStatisticsStartTime", 0);
    bool timeStatisticsOn = mesh.solutionDict().subDict("options").lookupOrDefault<bool>("timeStatisticsOn", false);
    scalar avgStartTime = mesh.solutionDict().subDict("options").lookupOrDefault<scalar>("avgStartTime", 0);
    scalar corrStartTime = mesh.solutionDict().subDict("options").lookupOrDefault<scalar>("corrStartTime", 0);

    // Read relaxation factor on the pressure gradient control.
    scalar alpha(mesh.solutionDict().subDict("options").lookupOrDefault<scalar>("alpha",1.0));

    // Read information about flow direction/wall-to-wall direction.
    word crossChannelDirection = mesh.solutionDict().subDict("options").lookupOrDefault<word>("crossChannelDirection", "y");
    word flowDirection = mesh.solutionDict().subDict("options").lookupOrDefault<word>("flowDirection", "x");

    // Convert the flow direction/wall-to-wall direction information to integers and check to
    // see if they are different.  If not, give a warning.
    int upIndex = 0;
    if (crossChannelDirection == "x")
    {
        upIndex = 0;
    }
    else if (crossChannelDirection == "y")
    {
        upIndex = 1;
    }
    else if (crossChannelDirection == "z")
    {
        upIndex = 2;
    }

    vector flowDirectionVector(vector::zero);
    if (flowDirection == "x")
    {
        flowDirectionVector.x() = 1.0;;
    }
    else if (flowDirection == "y")
    {
        flowDirectionVector.y() = 1.0;;
    }
    else if (flowDirection == "z")
    {
        flowDirectionVector.z() = 1.0;;
    }

    if (crossChannelDirection == flowDirection)
    {
        WarningIn("DNSPisoChannelFoam")
            << "The flow direction is set to " << endl
            << "the same as the cross channel " << endl
            << "averaging direction." << endl;
    }
