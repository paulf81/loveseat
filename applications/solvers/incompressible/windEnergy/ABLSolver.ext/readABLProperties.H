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

    // PROPERTIES CONCERNING THE SOURCE TERMS.

       // Specify the type of source to use for momentum and temperature.  The
       // possible types are "given" and "computed".  
       // - The "given" type means that the source values are directly given
       //   and the momentum and temperature fields will react accordingly.  
       // - The "computed" type means that the mean velocity and temperature
       //   are given and the source terms that maintain them are computed. 
       word momentumSourceType(ABLProperties.lookup("momentumSourceType"));
       word temperatureSourceType(ABLProperties.lookup("temperatureSourceType"));

       // The desired sources are given as a table of values versus height and
       // time that need not match with the mesh heights or simulation time steps.
       // Interpolation to the mesh heights and simulation times is done, so
       // specify if you want linearly or spline interpolation.
       word sourceInterpType(ABLProperties.lookup("sourceInterpType"));

       // Read in the heights at which the sources are given.
       List<scalar> sourceHeightsSpecified(ABLProperties.lookup("sourceHeights"));

       // Read in the table of momentum source vs. time and height.
       List<List<scalar> > sourceTableMomentumX(ABLProperties.lookup("sourceTableMomentumX"));
       List<List<scalar> > sourceTableMomentumY(ABLProperties.lookup("sourceTableMomentumY"));
       List<List<scalar> > sourceTableMomentumZ(ABLProperties.lookup("sourceTableMomentumZ"));

       // Read in the table of temperature source vs. time and height.
       List<List<scalar> > sourceTableTemperature(ABLProperties.lookup("sourceTableTemperature"));

       // Check sizes of source tables.
       bool sizeFlag = false;
       label nSourceHeights = sourceHeightsSpecified.size();
       label nSourceMomentumXTimes = sourceTableMomentumX.size();
       label nSourceMomentumYTimes = sourceTableMomentumY.size();
       label nSourceMomentumZTimes = sourceTableMomentumZ.size();
       label nSourceTemperatureTimes = sourceTableTemperature.size();
       forAll(sourceTableMomentumX,i)
       {
           if (sourceTableMomentumX[i].size()-1 != nSourceHeights)
           {
               sizeFlag = true;
           }
       }
       forAll(sourceTableMomentumY,i)
       {
           if (sourceTableMomentumY[i].size()-1 != nSourceHeights)
           {
               sizeFlag = true;
           }
       }
       forAll(sourceTableMomentumZ,i)
       {
           if (sourceTableMomentumZ[i].size()-1 != nSourceHeights)
           {
               sizeFlag = true;
           }
       }
       forAll(sourceTableTemperature,i)
       {
           if (sourceTableTemperature[i].size()-1 != nSourceHeights)
           {
               sizeFlag = true;
           }
       }


       // Break the source tables into interpolation tables.
       List<scalar> sourceMomentumXTimesSpecified(nSourceMomentumXTimes,0.0);
       List<scalar> sourceMomentumYTimesSpecified(nSourceMomentumYTimes,0.0);
       List<scalar> sourceMomentumZTimesSpecified(nSourceMomentumZTimes,0.0);
       List<scalar> sourceTemperatureTimesSpecified(nSourceTemperatureTimes,0.0);

       List<List<scalar> > sourceMomentumXSpecified(nSourceMomentumXTimes,List<scalar>(nSourceHeights,0.0));
       List<List<scalar> > sourceMomentumYSpecified(nSourceMomentumYTimes,List<scalar>(nSourceHeights,0.0));
       List<List<scalar> > sourceMomentumZSpecified(nSourceMomentumZTimes,List<scalar>(nSourceHeights,0.0));
       List<List<scalar> > sourceTemperatureSpecified(nSourceTemperatureTimes,List<scalar>(nSourceHeights,0.0));

       for(int i = 0, i < nSourceMomentumXTimes, i++)
       {
           sourceTimesSpecified[i] = sourceTableMomentumX[i][0];
           for(int j = 0, j < nH, j++)
           {
               sourceMomentumXSpecified[i][j] = sourceTableMomentumX[i][j+1];
               sourceMomentumYSpecified[i][j] = sourceTableMomentumY[i][j+1];
               sourceMomentumZSpecified[i][j] = sourceTableMomentumZ[i][j+1];
               sourceTemperatureSpecified[i][j] = sourceTableTemperature[i][j+1];
           }
       }










       // Read in a table of wind speed and direction versus time.
       List<List<scalar> > desiredWindTable(ABLProperties.lookup("desiredWindTable"));
       scalarField desiredWindTime(desiredWindTable.size(),0.0);
       scalarField desiredWindSpeed(desiredWindTable.size(),0.0);
       scalarField desiredWindDirectionDeg(desiredWindTable.size(),0.0);
       scalarField desiredWindDirection(desiredWindTable.size(),0.0);
       vectorField desiredWindVector(desiredWindTable.size(),vector::zero);
       forAll(desiredWindTime,i)
       {
           desiredWindTime[i]      = desiredWindTable[i][0];
           desiredWindSpeed[i]     = desiredWindTable[i][1];
           desiredWindDirectionDeg[i] = desiredWindTable[i][2];
       }

       // Height at which horizontally-averaged wind vector is specified
       dimensionedScalar hWind(ABLProperties.lookup("hWind"));

       // Write out the desired wind table to the log file.
       Info << "     Specified wind at " << hWind.value() << " m:" << endl;
       forAll(desiredWindTime,i)
       {
           Info << "              Time " << desiredWindTime[i] << "," << tab << "wind from " << desiredWindDirectionDeg[i] << " degrees at " << desiredWindSpeed[i] << " m/s" << endl;
       }
       Info << endl;
       Info << "     +x is east and +y is north" << endl;
       Info << "                   N" << endl;
       Info << "                   0" << endl;
       Info << "                   |" << endl << endl;
       Info << "        W 270 --       --  90 E" << endl << endl;
       Info << "                   |" << endl;
       Info << "                  180" << endl;
       Info << "                   S"  << endl;
      
       // Convert the cardinal wind direction (like on a compass) to the normal
       // convention of 0 radians on the + x axis with positive degrees in the
       // counter-clockwise direction.  Also, calculate the wind vector.
       forAll(desiredWindTime,i)
       {
           if (desiredWindDirectionDeg[i] > 180.0)
           {
               desiredWindDirection[i] = desiredWindDirectionDeg[i] - 180.0;
           }
           else
           {
               desiredWindDirection[i] = desiredWindDirectionDeg[i] + 180.0;
           }
           desiredWindDirection[i] = 90.0 - desiredWindDirection[i];
           if (desiredWindDirection[i] < 0.0)
           {
               desiredWindDirection[i] = desiredWindDirection[i] + 360.0;
           }
           desiredWindDirection[i] = desiredWindDirection[i] * ((Foam::constant::mathematical::pi)/180.0);

           vector windSpeedToVector;
           windSpeedToVector.x() = Foam::cos(desiredWindDirection[i]);
           windSpeedToVector.y() = Foam::sin(desiredWindDirection[i]);
           windSpeedToVector.z() = 0.0;
           desiredWindVector[i] = desiredWindSpeed[i] * windSpeedToVector;
       }

       // Relaxation factor on the pressure gradient control
       scalar alpha(ABLProperties.lookupOrDefault<scalar>("alpha",1.0));




    // PROPERTIES CONCERNING CORIOLIS FORCES

       // Planetary rotation period (hours)
       scalar planetaryRotationPeriod(readScalar(ABLProperties.lookup("planetaryRotationPeriod")));

       // Latitude on the planetary body (degrees)
       scalar latitude(readScalar(ABLProperties.lookup("latitude")));

       // Up index
       label upIndex = 2;
       vector nUp(vector::zero);
       nUp.z() = 1.0;

       // Compute the planetar rotation vector
       vector Omega_;
       Omega_.x() = 0.0;
       Omega_.y() = ((2.0 * Foam::constant::mathematical::pi) / (max(1.0E-5,planetaryRotationPeriod)*3600.0)) * Foam::cos(latitude*Foam::constant::mathematical::pi/180.0);
       Omega_.z() = ((2.0 * Foam::constant::mathematical::pi) / (max(1.0E-5,planetaryRotationPeriod)*3600.0)) * Foam::sin(latitude*Foam::constant::mathematical::pi/180.0);
       uniformDimensionedVectorField Omega
       (
           IOobject
           (
               "Omega",
               runTime.constant(),
               mesh,
               IOobject::NO_READ,
               IOobject::NO_WRITE
           ),
           dimensionedVector("Omega",dimensionSet(0, 0, -1, 0, 0, 0, 0),Omega_)
       );

       Info << Omega << endl;       




    // PROPERTIES CONCERNING GATHERING STATISTICS

       // Gather/write statistics?
       bool statisticsOn(ABLProperties.lookupOrDefault<bool>("statisticsOn", false));

       // Statistics gathering/writing frequency?
       int statisticsFreq(int(readScalar(ABLProperties.lookup("statisticsFrequency"))));




    // Name of the bottom boundary
    word lowerBoundaryName(ABLProperties.lookupOrDefault<word>("lowerBoundaryName","lower"));

    // Name of the upper boundary
    word upperBoundaryName(ABLProperties.lookupOrDefault<word>("upperBoundaryName","upper"));
