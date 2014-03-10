    IOdictionary ABLProperties
    (
        IOobject
        (
            "ABLProperties",
            runTime.time().constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Drive wind to specified speed a specified height using time
    // varying driving pressure gradient?
    bool driveWindOn(ABLProperties.lookupOrDefault<bool>("driveWindOn", false));

    // Desired horizontally-averaged wind speed at a height
    dimensionedScalar UWindSpeed(ABLProperties.lookup("UWindSpeed"));

    // Desired horizontally-averaged wind direction at a height (degrees)
    scalar UWindDir(ABLProperties.lookupOrDefault<scalar>("UWindDir",270.0));

    // Height at which horizontally-averaged wind vector is specified
    dimensionedScalar hWind(ABLProperties.lookup("hWind"));

    // Convert the cardinal wind direction (like on a compass) to the normal
    // convention of 0 radians on the + x axis with positive degrees in the
    // counter-clockwise direction
    Info << "     Specified wind at " << hWind.value() << " m is from " << UWindDir << " degrees at " << UWindSpeed.value() << " m/s" << endl;
    if (UWindDir > 180.0)
    {
        UWindDir = UWindDir - 180.0;
    }
    else
    {
        UWindDir = UWindDir + 180.0;
    }
    UWindDir = 90.0 - UWindDir;
    if (UWindDir < 0.0)
    {
        UWindDir = UWindDir + 360.0;
    }
    UWindDir = UWindDir * ((Foam::constant::mathematical::pi)/180.0);

    // Calculate the wind vector
    vector UWindSpeedToVector;
    UWindSpeedToVector.x() = Foam::cos(UWindDir);
    UWindSpeedToVector.y() = Foam::sin(UWindDir);
    UWindSpeedToVector.z() = 0.0;
    dimensionedVector UWind = UWindSpeed * UWindSpeedToVector;
    Info << "     This is a wind vector of " << UWind.value() << " m/s, where +x is east and +y is north" << endl;
    Info << "                               N" << endl;
    Info << "                               0" << endl;
    Info << "                               |" << endl << endl;
    Info << "                    W 270 --       --  90 E" << endl << endl;
    Info << "                               |" << endl;
    Info << "                              180" << endl;
    Info << "                               S"  << endl;


    // Relaxation factor on the pressure gradient control
    scalar alpha(ABLProperties.lookupOrDefault<scalar>("alpha",1.0));

    // Name of the bottom boundary
    word lowerBoundaryName(ABLProperties.lookupOrDefault<word>("lowerBoundaryName","lower"));

    // Name of the upper boundary
    word upperBoundaryName(ABLProperties.lookupOrDefault<word>("upperBoundaryName","upper"));

    // Up index
    label upIndex = 2;
    vector nUp(vector::zero);
    nUp.z() = 1.0;

    // Gather/write statistics?
    bool statisticsOn(ABLProperties.lookupOrDefault<bool>("statisticsOn", false));

    // Statistics gathering/writing frequency?
    int statisticsFreq(ABLProperties.lookupOrDefault<int>("statisticsFrequency", 0));

