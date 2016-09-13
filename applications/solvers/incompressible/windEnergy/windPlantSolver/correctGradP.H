  {
        volScalarField rAU("rAU", 1.0/UEqn.A());
        surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

        // Calculate the average wind velocity in the height range
        vector UWindI = vector::zero;
        forAll (windHeightCell, celli)
        {
            UWindI += U[windHeightCell[celli]] * mesh.V()[windHeightCell[celli]];
        }
        reduce(UWindI,sumOp<vector>());
        UWindI = UWindI / windHeightCellVolumeTotal;

        dimensionedVector UWindStar
	(
	    "UWindStar",
	    dimensionSet(0, 1, -1, 0, 0, 0, 0),
	    UWindI
	);

        dimensionedVector UWindStarParallel = UWindStar - ((UWindStar & nUp) * nUp);

        // Use the look-up table of desired wind speed to get the current wind vector.
        vector UWind_(vector::zero);
        scalar UWindSpeed_(0.0);
        scalar UWindDir_(0.0);

        if (windForcingType == "givenWindVector")
        {
            if (forcingInterpType == "linear")
            {
                UWind_      = interpolateXY(runTime.time().value(),desiredWindTime,desiredWindVector);
                UWindSpeed_ = interpolateXY(runTime.time().value(),desiredWindTime,desiredWindSpeed);
                UWindDir_   = interpolateXY(runTime.time().value(),desiredWindTime,desiredWindDirectionDeg);
            }
            else if (forcingInterpType == "cubic")
            {
                UWind_      = interpolateSplineXY(runTime.time().value(),desiredWindTime,desiredWindVector);
                UWindSpeed_ = interpolateSplineXY(runTime.time().value(),desiredWindTime,desiredWindSpeed);
                UWindDir_   = interpolateSplineXY(runTime.time().value(),desiredWindTime,desiredWindDirectionDeg);
            }
        }

        dimensionedVector UWind
        (
            "UWind",
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            UWind_
        );


        // Use the look-up table of desired driving pressure gradient to get the current value.
        vector gradPInterp(vector::zero);
        if (windForcingType == "givenPressureGradient")
        {
            if (forcingInterpType == "linear")
            {
                gradPInterp = interpolateXY(runTime.time().value(),desiredGradPTime,desiredGradPVector);
            }
            else if (forcingInterpType == "cubic")
            {
                gradPInterp = interpolateSplineXY(runTime.time().value(),desiredGradPTime,desiredGradPVector);
            }
        }
        
        dimensionedVector gradP_
        (
            "gradP_",
            dimensionSet(0, 1, -2, 0, 0, 0, 0),
            gradPInterp
        );
        gradP = gradP_;


       
        // Correct the driving pressure gradient and wind 
        if (windForcingType == "givenWindVector")
        {
 	     dimensionedVector pDelta = alpha*(UWind - UWindStarParallel)/rAU.weightedAverage(mesh.V());

 	     U += rAU*pDelta;
	     phi += rAUf*(pDelta & mesh.Sf());

	     gradP -=  ((gradP & nUp) * nUp) + pDelta;
        }
        else if (windForcingType == "givenPressureGradient")
        {
        }

	Info << "Uncorrected UWind = "  << UWindStar.value() << tab
             << "Pressure Gradient = " << gradP.value()   << endl;
  }
