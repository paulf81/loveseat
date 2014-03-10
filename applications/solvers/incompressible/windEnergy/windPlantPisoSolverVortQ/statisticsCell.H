   if (statisticsOn)
   {
        if (timeStep % statisticsFreq == 0)
	{
	     // Average the field variables
             #include "averageFields.H"

	     // Then get the statistics at each vertical level
	     forAll(hLevelsValues,hLevelsI)
	     {
	          symmTensor velFluxMeanVol = symmTensor::zero;
		  symmTensor velFluxFluxMeanVol = symmTensor::zero;
		  vector tempFluxMeanVol = vector::zero;

	  	  // sum variable times cell volume over the cells at a certain level
		  for (label i = 0; i < numCellPerLevel[hLevelsI]; i++)
		  {
		       label cellI = hLevelsCellList[hLevelsI][i];
	               velFluxMeanVol.xx() += Uprime[cellI].x() * Uprime[cellI].x() * mesh.V()[cellI];
	               velFluxMeanVol.xy() += Uprime[cellI].x() * Uprime[cellI].y() * mesh.V()[cellI];
	               velFluxMeanVol.xz() += Uprime[cellI].x() * Uprime[cellI].z() * mesh.V()[cellI];
	               velFluxMeanVol.yy() += Uprime[cellI].y() * Uprime[cellI].y() * mesh.V()[cellI];
	               velFluxMeanVol.yz() += Uprime[cellI].y() * Uprime[cellI].z() * mesh.V()[cellI];
	               velFluxMeanVol.zz() += Uprime[cellI].z() * Uprime[cellI].z() * mesh.V()[cellI];

	               velFluxFluxMeanVol.xx() += Uprime[cellI].z() * Uprime[cellI].x() * Uprime[cellI].x() * mesh.V()[cellI];
	               velFluxFluxMeanVol.xy() += Uprime[cellI].z() * Uprime[cellI].x() * Uprime[cellI].y() * mesh.V()[cellI];
	               velFluxFluxMeanVol.xz() += Uprime[cellI].z() * Uprime[cellI].x() * Uprime[cellI].z() * mesh.V()[cellI];
	               velFluxFluxMeanVol.yy() += Uprime[cellI].z() * Uprime[cellI].y() * Uprime[cellI].y() * mesh.V()[cellI];
	               velFluxFluxMeanVol.yz() += Uprime[cellI].z() * Uprime[cellI].y() * Uprime[cellI].z() * mesh.V()[cellI];
	               velFluxFluxMeanVol.zz() += Uprime[cellI].z() * Uprime[cellI].z() * Uprime[cellI].z() * mesh.V()[cellI];

		       tempFluxMeanVol.x() += Tprime[cellI]     * Uprime[cellI].x() * mesh.V()[cellI];
		       tempFluxMeanVol.y() += Tprime[cellI]     * Uprime[cellI].y() * mesh.V()[cellI];
		       tempFluxMeanVol.z() += Tprime[cellI]     * Uprime[cellI].z() * mesh.V()[cellI];
		  }

		  // parallel gather the sums from each processor and scatter back out to processors
		  reduce(velFluxMeanVol,sumOp<symmTensor>());
		  reduce(velFluxFluxMeanVol,sumOp<symmTensor>());
		  reduce(tempFluxMeanVol,sumOp<vector>());

		  // divide the volume-weighted sums by total volume at a certain level
                  velFluxLevelsList[hLevelsI] = velFluxMeanVol/totVolPerLevel[hLevelsI];
                  velFluxFluxLevelsList[hLevelsI] = velFluxFluxMeanVol/totVolPerLevel[hLevelsI];
                  tempFluxLevelsList[hLevelsI] = tempFluxMeanVol/totVolPerLevel[hLevelsI];
	     }

             // Write the statistics to files
	     if (Pstream::myProcNo() == 0)
             {
	          TmeanFile << runTime.timeName() << " " << runTime.deltaT().value();

	          umeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          vmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          wmeanFile << runTime.timeName() << " " << runTime.deltaT().value();

	          uuPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          uvPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          uwPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          vvPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          vwPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          wwPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();

	          wuuPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          wuvPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          wuwPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          wvvPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          wvwPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          wwwPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();

	          TuPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          TvPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          TwPmeanFile << runTime.timeName() << " " << runTime.deltaT().value();

	          forAll(hLevelsValues,hLevelsI)
	          {
	               TmeanFile << " " << TmeanLevelsList[hLevelsI];

	               umeanFile << " " << UmeanLevelsList[hLevelsI].x();
	               vmeanFile << " " << UmeanLevelsList[hLevelsI].y();
	               wmeanFile << " " << UmeanLevelsList[hLevelsI].z();

	               uuPmeanFile << " " << velFluxLevelsList[hLevelsI].xx();
	               uvPmeanFile << " " << velFluxLevelsList[hLevelsI].xy();
	               uwPmeanFile << " " << velFluxLevelsList[hLevelsI].xz();
	               vvPmeanFile << " " << velFluxLevelsList[hLevelsI].yy();
	               vwPmeanFile << " " << velFluxLevelsList[hLevelsI].yz();
	               wwPmeanFile << " " << velFluxLevelsList[hLevelsI].zz();

	               wuuPmeanFile << " " << velFluxFluxLevelsList[hLevelsI].xx();
	               wuvPmeanFile << " " << velFluxFluxLevelsList[hLevelsI].xy();
	               wuwPmeanFile << " " << velFluxFluxLevelsList[hLevelsI].xz();
	               wvvPmeanFile << " " << velFluxFluxLevelsList[hLevelsI].yy();
	               wvwPmeanFile << " " << velFluxFluxLevelsList[hLevelsI].yz();
	               wwwPmeanFile << " " << velFluxFluxLevelsList[hLevelsI].zz();

	               TuPmeanFile << " " << tempFluxLevelsList[hLevelsI].x();
	               TvPmeanFile << " " << tempFluxLevelsList[hLevelsI].y();
	               TwPmeanFile << " " << tempFluxLevelsList[hLevelsI].z();
	          }

	          TmeanFile << endl;

	          umeanFile << endl;
	          vmeanFile << endl;
	          wmeanFile << endl;

	          uuPmeanFile << endl;
	          uvPmeanFile << endl;
	          uwPmeanFile << endl;
	          vvPmeanFile << endl;
	          vwPmeanFile << endl;
	          wwPmeanFile << endl;

	          wuuPmeanFile << endl;
	          wuvPmeanFile << endl;
	          wuwPmeanFile << endl;
	          wvvPmeanFile << endl;
	          wvwPmeanFile << endl;
	          wwwPmeanFile << endl;

	          TuPmeanFile << endl;
	          TvPmeanFile << endl;
	          TwPmeanFile << endl;
	     }
        }
   }
