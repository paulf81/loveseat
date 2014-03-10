  {
        // Calculate the average wind velocity at the closest level to specified
        //vector UWind1 = vector::zero;
        //for (label i = 0; i < numCellPerLevel[hLevelsWind1I]; i++)
        //{
	//     label cellI = hLevelsCellList[hLevelsWind1I][i];
	//     UWind1 += U[cellI] * mesh.V()[cellI];
        //}
        //reduce(UWind1,sumOp<vector>());
	//UWind1 = UWind1/totVolPerLevel[hLevelsWind1I];

        // Calculate the average wind velocity at the next closest level to specified
	//vector UWind2 = vector::zero;
        //for (label i = 0; i < numCellPerLevel[hLevelsWind2I]; i++)
        //{
	//     label cellI = hLevelsCellList[hLevelsWind2I][i];
	//     UWind2 += U[cellI] * mesh.V()[cellI];
        //}
        //reduce(UWind2,sumOp<vector>());
	//UWind2 = UWind2/totVolPerLevel[hLevelsWind2I];


        //dimensionedVector UWindStar
	//(
	//    "UWindStar",
	//    dimensionSet(0, 1, -1, 0, 0, 0, 0),
	//    UWind1 + (((UWind2 - UWind1)/(hLevelsWind2 - hLevelsWind1)) * (hWind.value() - hLevelsWind1))
	//);
        
        //if (driveWindOn)
        //{
        //     dimensionedVector UWindStarParallel = UWindStar - ((UWindStar & nUp) * nUp);

	//     dimensionedVector pdDelta = alpha*(UWind - UWindStarParallel)/rUA.weightedAverage(mesh.V());

 	//     U += rUA*pdDelta;
	//     phi += rUAf*(pdDelta & mesh.Sf());

	//     gradPd -=  ((gradPd & nUp) * nUp) + pdDelta;
        //}

	//Info << "Uncorrected UWind = "  << UWindStar.value() << tab
        //     << "Pressure Gradient = " << gradPd.value()   << endl;
        Info << "Pressure Gradient = " << gradPd.value() << endl;
  }
