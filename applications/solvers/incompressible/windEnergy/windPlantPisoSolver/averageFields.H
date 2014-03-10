    // Get horizontally averaged velocity profile
    //forAll(hLevelsValues,hLevelsI)
    //{
    //    scalar TmeanVol = 0.0;
    //    vector UmeanVol = vector::zero;
    //    for (label i = 0; i < numCellPerLevel[hLevelsI]; i++)
    //    {
//	    label cellI = hLevelsCellList[hLevelsI][i];
//	    TmeanVol += T[cellI] * mesh.V()[cellI];
//	    UmeanVol += U[cellI] * mesh.V()[cellI];
    //    }
    //    reduce(TmeanVol,sumOp<scalar>());
    //    reduce(UmeanVol,sumOp<vector>());
    //    TmeanLevelsList[hLevelsI] = TmeanVol/totVolPerLevel[hLevelsI];
    //    UmeanLevelsList[hLevelsI] = UmeanVol/totVolPerLevel[hLevelsI];


    //    for(label i = 0; i < numCellPerLevel[hLevelsI]; i++)
//	{
//	    label cellI = hLevelsCellList[hLevelsI][i];
//	    Tmean[cellI] = TmeanLevelsList[hLevelsI];
//	    Umean[cellI] = UmeanLevelsList[hLevelsI];
//	}
    //}


   
    if (runTime.value() >= meanAvgStartTime)
    {
        Umean = (Umean * meanTimeSum) + (runTime.deltaT().value() * U);
        Tmean = (Tmean * meanTimeSum) + (runTime.deltaT().value() * T);

        meanTimeSum += runTime.deltaT().value();

        Umean = Umean / meanTimeSum;
        Tmean = Tmean / meanTimeSum;

        Uprime = U - Umean;
        Tprime = T - Tmean;
    }

    if (runTime.value() >= corrAvgStartTime)
    {
        UprimeUprimeMean = (UprimeUprimeMean * corrTimeSum) + (runTime.deltaT().value() * sqr(Uprime));

        TprimeUprimeMean  = (TprimeUprimeMean  * corrTimeSum) + (runTime.deltaT().value() * Tprime * Uprime);

        corrTimeSum += runTime.deltaT().value();

        UprimeUprimeMean = UprimeUprimeMean / corrTimeSum;
        TprimeUprimeMean = TprimeUprimeMean / corrTimeSum;
    }

    k = 0.5 * tr(UprimeUprimeMean);
