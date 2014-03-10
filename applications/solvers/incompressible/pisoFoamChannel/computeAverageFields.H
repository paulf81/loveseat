{
    volVectorField uPrime = U;
    volScalarField pPrime = p;

    // Compute the running time averages of the fields and then get the fluctuation.
    if (runTime.value() >= avgStartTime)
    {
        UAvg = (UAvg * avgTimeSum) + (runTime.deltaT().value() * U);
        pAvg = (pAvg * avgTimeSum) + (runTime.deltaT().value() * p);

        avgTimeSum += runTime.deltaT().value();

        UAvg = UAvg / avgTimeSum;
        pAvg = pAvg / avgTimeSum;
    }

    // Compute the running time averages of the correlation fields.
    if (runTime.value() >= corrStartTime)
    {
        uPrime = U - UAvg;
        pPrime = p - pAvg;

        uuPrime2 = (uuPrime2 * corrTimeSum) + (runTime.deltaT().value() * Foam::sqr(uPrime));
        upPrime2 = (upPrime2 * corrTimeSum) + (runTime.deltaT().value() * pPrime * uPrime);
        ppPrime2 = (ppPrime2 * corrTimeSum) + (runTime.deltaT().value() * Foam::sqr(pPrime));

        uuuPrime3 = (uuuPrime3 * corrTimeSum) + (runTime.deltaT().value() * Foam::sqr(uPrime) * uPrime.component(0));
        vuuPrime3 = (vuuPrime3 * corrTimeSum) + (runTime.deltaT().value() * Foam::sqr(uPrime) * uPrime.component(1));
        wuuPrime3 = (wuuPrime3 * corrTimeSum) + (runTime.deltaT().value() * Foam::sqr(uPrime) * uPrime.component(2));

        corrTimeSum += runTime.deltaT().value();

        uuPrime2 = uuPrime2 / corrTimeSum;
        upPrime2 = upPrime2 / corrTimeSum;
        ppPrime2 = ppPrime2 / corrTimeSum;
        uuuPrime3 = uuuPrime3 / corrTimeSum;
        vuuPrime3 = vuuPrime3 / corrTimeSum;
        wuuPrime3 = wuuPrime3 / corrTimeSum;

        forAll(uRMS,cellI)
        {
            uRMS[cellI].x() = Foam::sqrt(uuPrime2[cellI].xx());
            uRMS[cellI].y() = Foam::sqrt(uuPrime2[cellI].yy());
            uRMS[cellI].z() = Foam::sqrt(uuPrime2[cellI].zz());
        }
        pRMS = Foam::sqrt(ppPrime2);
    }

    // If at the end of the simulation, horizontally average the time averages
    // to get a profile.
    if (runTime.run() == false)
    {
        // Define a set of mean variable lists.
        List<vector>     UAvgLevelsList(hLevelsTotal,vector::zero);
        List<scalar>     pAvgLevelsList(hLevelsTotal,0.0);

        List<vector>     uRMSLevelsList(hLevelsTotal,vector::zero);
        List<scalar>     pRMSLevelsList(hLevelsTotal,0.0);

        List<symmTensor> uuPrime2LevelsList(hLevelsTotal,symmTensor::zero);
        List<vector>     upPrime2LevelsList(hLevelsTotal,vector::zero);

        List<symmTensor> uuuPrime3LevelsList(hLevelsTotal,symmTensor::zero);
        List<symmTensor> vuuPrime3LevelsList(hLevelsTotal,symmTensor::zero);
        List<symmTensor> wuuPrime3LevelsList(hLevelsTotal,symmTensor::zero);

        // Get the average profiles.
        forAll(hLevelsValues,hLevelsI)
        {
            vector UAvgVol = vector::zero;
            scalar pAvgVol = 0.0;

            vector uRMSVol = vector::zero;
            scalar pRMSVol = 0.0;

            symmTensor uuPrime2Vol = symmTensor::zero;
            vector upPrime2Vol = vector::zero;
            scalar ppPrime2Vol = 0.0;

            symmTensor uuuPrime3Vol = symmTensor::zero;
            symmTensor vuuPrime3Vol = symmTensor::zero;
            symmTensor wuuPrime3Vol = symmTensor::zero;

            // sum variable times cell volume over the cells at a certain level.
            for (label i = 0; i < numCellPerLevel[hLevelsI]; i++)
            {
                label cellI = hLevelsCellList[hLevelsI][i];
                
                UAvgVol += UAvg[cellI] * mesh.V()[cellI];
                pAvgVol += pAvg[cellI] * mesh.V()[cellI];

                uRMSVol += uRMS[cellI] * mesh.V()[cellI];
                pRMSVol += pRMS[cellI] * mesh.V()[cellI];

                uuPrime2Vol += uuPrime2[cellI] * mesh.V()[cellI];
                upPrime2Vol += upPrime2[cellI] * mesh.V()[cellI];
                ppPrime2Vol += ppPrime2[cellI] * mesh.V()[cellI];

                uuuPrime3Vol += uuuPrime3[cellI] * mesh.V()[cellI];
                vuuPrime3Vol += vuuPrime3[cellI] * mesh.V()[cellI];
                wuuPrime3Vol += wuuPrime3[cellI] * mesh.V()[cellI];
            }
            reduce(UAvgVol,sumOp<vector>());
            reduce(pAvgVol,sumOp<scalar>());
            reduce(uRMSVol,sumOp<vector>());
            reduce(pRMSVol,sumOp<scalar>());
            reduce(uuPrime2Vol,sumOp<symmTensor>());
            reduce(upPrime2Vol,sumOp<vector>());
            reduce(ppPrime2Vol,sumOp<scalar>());
            reduce(uuuPrime3Vol,sumOp<symmTensor>());
            reduce(vuuPrime3Vol,sumOp<symmTensor>());
            reduce(wuuPrime3Vol,sumOp<symmTensor>());

            UAvgLevelsList[hLevelsI] = UAvgVol/totVolPerLevel[hLevelsI];
            pAvgLevelsList[hLevelsI] = pAvgVol/totVolPerLevel[hLevelsI];
            uRMSLevelsList[hLevelsI] = uRMSVol/totVolPerLevel[hLevelsI];
            pRMSLevelsList[hLevelsI] = pRMSVol/totVolPerLevel[hLevelsI];
            uuPrime2LevelsList[hLevelsI] = uuPrime2Vol/totVolPerLevel[hLevelsI];
            upPrime2LevelsList[hLevelsI] = upPrime2Vol/totVolPerLevel[hLevelsI];
            uuuPrime3LevelsList[hLevelsI] = uuuPrime3Vol/totVolPerLevel[hLevelsI];
            vuuPrime3LevelsList[hLevelsI] = vuuPrime3Vol/totVolPerLevel[hLevelsI];
            wuuPrime3LevelsList[hLevelsI] = wuuPrime3Vol/totVolPerLevel[hLevelsI];
        }


        // Write the statistics to files.
        if (Pstream::myProcNo() == 0)
        {
             UMeanFile2 << runTime.timeName();
             VMeanFile2 << runTime.timeName();
             WMeanFile2 << runTime.timeName();
             pMeanFile2 << runTime.timeName();

             uRMSFile2 << runTime.timeName();
             vRMSFile2 << runTime.timeName();
             wRMSFile2 << runTime.timeName();
             pRMSFile2 << runTime.timeName();

             uuPrime2File2 << runTime.timeName();
             uvPrime2File2 << runTime.timeName();
             uwPrime2File2 << runTime.timeName();
             vvPrime2File2 << runTime.timeName();
             vwPrime2File2 << runTime.timeName();
             wwPrime2File2 << runTime.timeName();

             upPrime2File2 << runTime.timeName();
             vpPrime2File2 << runTime.timeName();
             wpPrime2File2 << runTime.timeName();

             uuuPrime3File2 << runTime.timeName();
             uuvPrime3File2 << runTime.timeName();
             uuwPrime3File2 << runTime.timeName();
             uvvPrime3File2 << runTime.timeName();
             uvwPrime3File2 << runTime.timeName();
             uwwPrime3File2 << runTime.timeName();

             vuuPrime3File2 << runTime.timeName();
             vuvPrime3File2 << runTime.timeName();
             vuwPrime3File2 << runTime.timeName();
             vvvPrime3File2 << runTime.timeName();
             vvwPrime3File2 << runTime.timeName();
             vwwPrime3File2 << runTime.timeName();

             wuuPrime3File2 << runTime.timeName();
             wuvPrime3File2 << runTime.timeName();
             wuwPrime3File2 << runTime.timeName();
             wvvPrime3File2 << runTime.timeName();
             wvwPrime3File2 << runTime.timeName();
             wwwPrime3File2 << runTime.timeName();


             forAll(hLevelsValues,hLevelsI)
             {
                 UMeanFile2 << " " << UAvgLevelsList[hLevelsI].x();
                 VMeanFile2 << " " << UAvgLevelsList[hLevelsI].y();
                 WMeanFile2 << " " << UAvgLevelsList[hLevelsI].z();
                 pMeanFile2 << " " << pAvgLevelsList[hLevelsI];

                 uRMSFile2 << " " << uRMSLevelsList[hLevelsI].x();
                 vRMSFile2 << " " << uRMSLevelsList[hLevelsI].y();
                 wRMSFile2 << " " << uRMSLevelsList[hLevelsI].z();
                 pRMSFile2 << " " << pRMSLevelsList[hLevelsI];

                 uuPrime2File2 << " " << uuPrime2LevelsList[hLevelsI].xx();
                 uvPrime2File2 << " " << uuPrime2LevelsList[hLevelsI].xy();
                 uwPrime2File2 << " " << uuPrime2LevelsList[hLevelsI].xz();
                 vvPrime2File2 << " " << uuPrime2LevelsList[hLevelsI].yy();
                 vwPrime2File2 << " " << uuPrime2LevelsList[hLevelsI].yz();
                 wwPrime2File2 << " " << uuPrime2LevelsList[hLevelsI].zz();

                 upPrime2File2 << " " << upPrime2LevelsList[hLevelsI].x();
                 vpPrime2File2 << " " << upPrime2LevelsList[hLevelsI].y();
                 wpPrime2File2 << " " << upPrime2LevelsList[hLevelsI].z();

                 uuuPrime3File2 << " " << uuuPrime3LevelsList[hLevelsI].xx();
                 uuvPrime3File2 << " " << uuuPrime3LevelsList[hLevelsI].xy();
                 uuwPrime3File2 << " " << uuuPrime3LevelsList[hLevelsI].xz();
                 uvvPrime3File2 << " " << uuuPrime3LevelsList[hLevelsI].yy();
                 uvwPrime3File2 << " " << uuuPrime3LevelsList[hLevelsI].yz();
                 uwwPrime3File2 << " " << uuuPrime3LevelsList[hLevelsI].zz();

                 vuuPrime3File2 << " " << vuuPrime3LevelsList[hLevelsI].xx();
                 vuvPrime3File2 << " " << vuuPrime3LevelsList[hLevelsI].xy();
                 vuwPrime3File2 << " " << vuuPrime3LevelsList[hLevelsI].xz();
                 vvvPrime3File2 << " " << vuuPrime3LevelsList[hLevelsI].yy();
                 vvwPrime3File2 << " " << vuuPrime3LevelsList[hLevelsI].yz();
                 vwwPrime3File2 << " " << vuuPrime3LevelsList[hLevelsI].zz();

                 wuuPrime3File2 << " " << wuuPrime3LevelsList[hLevelsI].xx();
                 wuvPrime3File2 << " " << wuuPrime3LevelsList[hLevelsI].xy();
                 wuwPrime3File2 << " " << wuuPrime3LevelsList[hLevelsI].xz();
                 wvvPrime3File2 << " " << wuuPrime3LevelsList[hLevelsI].yy();
                 wvwPrime3File2 << " " << wuuPrime3LevelsList[hLevelsI].yz();
                 wwwPrime3File2 << " " << wuuPrime3LevelsList[hLevelsI].zz();
             }

             UMeanFile2 << endl;
             VMeanFile2 << endl;
             WMeanFile2 << endl;
             pMeanFile2 << endl;
             uRMSFile2 << endl;
             vRMSFile2 << endl;
             wRMSFile2 << endl;
             pRMSFile2 << endl;

             uuPrime2File2 << endl;
             uvPrime2File2 << endl;
             uwPrime2File2 << endl;
             vvPrime2File2 << endl;
             vwPrime2File2 << endl;
             wwPrime2File2 << endl;

             upPrime2File2 << endl;
             vpPrime2File2 << endl;
             wpPrime2File2 << endl;

             uuuPrime3File2 << endl;
             uuvPrime3File2 << endl;
             uuwPrime3File2 << endl;
             uvvPrime3File2 << endl;
             uvwPrime3File2 << endl;
             uwwPrime3File2 << endl;

             vuuPrime3File2 << endl;
             vuvPrime3File2 << endl;
             vuwPrime3File2 << endl;
             vvvPrime3File2 << endl;
             vvwPrime3File2 << endl;
             vwwPrime3File2 << endl;

             wuuPrime3File2 << endl;
             wuvPrime3File2 << endl;
             wuwPrime3File2 << endl;
             wvvPrime3File2 << endl;
             wvwPrime3File2 << endl;
             wwwPrime3File2 << endl;
        }   
    }
}
