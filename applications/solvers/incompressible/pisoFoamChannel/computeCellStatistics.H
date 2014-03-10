        {
             // Define a set of mean variable lists.
             List<vector>     UMeanLevelsList(hLevelsTotal,vector::zero);
             List<scalar>     pMeanLevelsList(hLevelsTotal,0.0);

             List<vector>     uRMSLevelsList(hLevelsTotal,vector::zero);
             List<scalar>     pRMSLevelsList(hLevelsTotal,0.0);

             List<symmTensor> uuPrime2LevelsList(hLevelsTotal,symmTensor::zero);
             List<vector>     upPrime2LevelsList(hLevelsTotal,vector::zero);

             List<symmTensor> uuuPrime3LevelsList(hLevelsTotal,symmTensor::zero);
             List<symmTensor> vuuPrime3LevelsList(hLevelsTotal,symmTensor::zero);
             List<symmTensor> wuuPrime3LevelsList(hLevelsTotal,symmTensor::zero);

             // Get horizontally averaged velocity profile.
             volVectorField UMean = U;
             volScalarField pMean = p;
             
             forAll(hLevelsValues,hLevelsI)
             {
                 vector UMeanVol = vector::zero;
                 scalar pMeanVol = 0.0;
                 for (label i = 0; i < numCellPerLevel[hLevelsI]; i++)
                 {
                     label cellI = hLevelsCellList[hLevelsI][i];
                     UMeanVol += U[cellI] * mesh.V()[cellI];
                     pMeanVol += p[cellI] * mesh.V()[cellI];
                 }
                 reduce(UMeanVol,sumOp<vector>());
                 reduce(pMeanVol,sumOp<scalar>());

                 UMeanLevelsList[hLevelsI] = UMeanVol/totVolPerLevel[hLevelsI];
                 pMeanLevelsList[hLevelsI] = pMeanVol/totVolPerLevel[hLevelsI];


                 for(label i = 0; i < numCellPerLevel[hLevelsI]; i++)
                 {
                     label cellI = hLevelsCellList[hLevelsI][i];
                     UMean[cellI] = UMeanLevelsList[hLevelsI];
                     pMean[cellI] = pMeanLevelsList[hLevelsI];
                 }
             }

             // Then get fluctuating part.
             volVectorField uPrime = U - UMean;
             volScalarField pPrime = p - pMean;

	     // Then get the statistics at each vertical level.
	     forAll(hLevelsValues,hLevelsI)
	     {
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
	               
                       uuPrime2Vol.xx() += uPrime[cellI].x() * uPrime[cellI].x() * mesh.V()[cellI];
	               uuPrime2Vol.xy() += uPrime[cellI].x() * uPrime[cellI].y() * mesh.V()[cellI];
	               uuPrime2Vol.xz() += uPrime[cellI].x() * uPrime[cellI].z() * mesh.V()[cellI];
	               uuPrime2Vol.yy() += uPrime[cellI].y() * uPrime[cellI].y() * mesh.V()[cellI];
	               uuPrime2Vol.yz() += uPrime[cellI].y() * uPrime[cellI].z() * mesh.V()[cellI];
	               uuPrime2Vol.zz() += uPrime[cellI].z() * uPrime[cellI].z() * mesh.V()[cellI];

	               upPrime2Vol.x() += uPrime[cellI].x() * pPrime[cellI] * mesh.V()[cellI];
	               upPrime2Vol.y() += uPrime[cellI].y() * pPrime[cellI] * mesh.V()[cellI];
	               upPrime2Vol.z() += uPrime[cellI].z() * pPrime[cellI] * mesh.V()[cellI];

	               ppPrime2Vol += pPrime[cellI] * pPrime[cellI] * mesh.V()[cellI];

	               uuuPrime3Vol.xx() += uPrime[cellI].x() * uPrime[cellI].x() * uPrime[cellI].x() * mesh.V()[cellI];
	               uuuPrime3Vol.xy() += uPrime[cellI].x() * uPrime[cellI].x() * uPrime[cellI].y() * mesh.V()[cellI];
	               uuuPrime3Vol.xz() += uPrime[cellI].x() * uPrime[cellI].x() * uPrime[cellI].z() * mesh.V()[cellI];
	               uuuPrime3Vol.yy() += uPrime[cellI].x() * uPrime[cellI].y() * uPrime[cellI].y() * mesh.V()[cellI];
	               uuuPrime3Vol.yz() += uPrime[cellI].x() * uPrime[cellI].y() * uPrime[cellI].z() * mesh.V()[cellI];
	               uuuPrime3Vol.zz() += uPrime[cellI].x() * uPrime[cellI].z() * uPrime[cellI].z() * mesh.V()[cellI];
	               
                       vuuPrime3Vol.xx() += uPrime[cellI].y() * uPrime[cellI].x() * uPrime[cellI].x() * mesh.V()[cellI];
	               vuuPrime3Vol.xy() += uPrime[cellI].y() * uPrime[cellI].x() * uPrime[cellI].y() * mesh.V()[cellI];
	               vuuPrime3Vol.xz() += uPrime[cellI].y() * uPrime[cellI].x() * uPrime[cellI].z() * mesh.V()[cellI];
	               vuuPrime3Vol.yy() += uPrime[cellI].y() * uPrime[cellI].y() * uPrime[cellI].y() * mesh.V()[cellI];
	               vuuPrime3Vol.yz() += uPrime[cellI].y() * uPrime[cellI].y() * uPrime[cellI].z() * mesh.V()[cellI];
	               vuuPrime3Vol.zz() += uPrime[cellI].y() * uPrime[cellI].z() * uPrime[cellI].z() * mesh.V()[cellI];
	               
                       wuuPrime3Vol.xx() += uPrime[cellI].z() * uPrime[cellI].x() * uPrime[cellI].x() * mesh.V()[cellI];
	               wuuPrime3Vol.xy() += uPrime[cellI].z() * uPrime[cellI].x() * uPrime[cellI].y() * mesh.V()[cellI];
	               wuuPrime3Vol.xz() += uPrime[cellI].z() * uPrime[cellI].x() * uPrime[cellI].z() * mesh.V()[cellI];
	               wuuPrime3Vol.yy() += uPrime[cellI].z() * uPrime[cellI].y() * uPrime[cellI].y() * mesh.V()[cellI];
	               wuuPrime3Vol.yz() += uPrime[cellI].z() * uPrime[cellI].y() * uPrime[cellI].z() * mesh.V()[cellI];
	               wuuPrime3Vol.zz() += uPrime[cellI].z() * uPrime[cellI].z() * uPrime[cellI].z() * mesh.V()[cellI];
		  }

		  // parallel gather the sums from each processor and scatter back out to processors.
		  reduce(uuPrime2Vol,sumOp<symmTensor>());
		  reduce(upPrime2Vol,sumOp<vector>());
		  reduce(ppPrime2Vol,sumOp<scalar>());
		  reduce(uuuPrime3Vol,sumOp<symmTensor>());
		  reduce(vuuPrime3Vol,sumOp<symmTensor>());
		  reduce(wuuPrime3Vol,sumOp<symmTensor>());

		  // divide the volume-weighted sums by total volume at a certain level.
                  uuPrime2LevelsList[hLevelsI] = uuPrime2Vol/totVolPerLevel[hLevelsI];
                  upPrime2LevelsList[hLevelsI] = upPrime2Vol/totVolPerLevel[hLevelsI];
                  uuuPrime3LevelsList[hLevelsI] = uuuPrime3Vol/totVolPerLevel[hLevelsI];
                  vuuPrime3LevelsList[hLevelsI] = vuuPrime3Vol/totVolPerLevel[hLevelsI];
                  wuuPrime3LevelsList[hLevelsI] = wuuPrime3Vol/totVolPerLevel[hLevelsI];

                  // Get the RMS values.
                  uRMSLevelsList[hLevelsI].x() = Foam::sqrt(uuPrime2LevelsList[hLevelsI].xx());
                  uRMSLevelsList[hLevelsI].y() = Foam::sqrt(uuPrime2LevelsList[hLevelsI].yy());
                  uRMSLevelsList[hLevelsI].z() = Foam::sqrt(uuPrime2LevelsList[hLevelsI].zz());
                  pRMSLevelsList[hLevelsI] = Foam::sqrt(ppPrime2Vol/totVolPerLevel[hLevelsI]);
	     }

             // Write the statistics to files.
	     if (Pstream::myProcNo() == 0)
             {
	          UMeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          VMeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          WMeanFile << runTime.timeName() << " " << runTime.deltaT().value();
	          pMeanFile << runTime.timeName() << " " << runTime.deltaT().value();

                  uRMSFile << runTime.timeName() << " " << runTime.deltaT().value();
                  vRMSFile << runTime.timeName() << " " << runTime.deltaT().value();
                  wRMSFile << runTime.timeName() << " " << runTime.deltaT().value();
                  pRMSFile << runTime.timeName() << " " << runTime.deltaT().value();

	          uuPrime2File << runTime.timeName() << " " << runTime.deltaT().value();
	          uvPrime2File << runTime.timeName() << " " << runTime.deltaT().value();
	          uwPrime2File << runTime.timeName() << " " << runTime.deltaT().value();
	          vvPrime2File << runTime.timeName() << " " << runTime.deltaT().value();
	          vwPrime2File << runTime.timeName() << " " << runTime.deltaT().value();
	          wwPrime2File << runTime.timeName() << " " << runTime.deltaT().value();
	          
                  upPrime2File << runTime.timeName() << " " << runTime.deltaT().value();
	          vpPrime2File << runTime.timeName() << " " << runTime.deltaT().value();
	          wpPrime2File << runTime.timeName() << " " << runTime.deltaT().value();

	          uuuPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          uuvPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          uuwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          uvvPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          uvwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          uwwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();

	          vuuPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          vuvPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          vuwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          vvvPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          vvwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          vwwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();

	          wuuPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          wuvPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          wuwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          wvvPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          wvwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();
	          wwwPrime3File << runTime.timeName() << " " << runTime.deltaT().value();


	          forAll(hLevelsValues,hLevelsI)
	          {
	               UMeanFile << " " << UMeanLevelsList[hLevelsI].x();
	               VMeanFile << " " << UMeanLevelsList[hLevelsI].y();
	               WMeanFile << " " << UMeanLevelsList[hLevelsI].z();
	               pMeanFile << " " << pMeanLevelsList[hLevelsI];

                       uRMSFile << " " << uRMSLevelsList[hLevelsI].x();
                       vRMSFile << " " << uRMSLevelsList[hLevelsI].y();
                       wRMSFile << " " << uRMSLevelsList[hLevelsI].z();
                       pRMSFile << " " << pRMSLevelsList[hLevelsI];

	               uuPrime2File << " " << uuPrime2LevelsList[hLevelsI].xx();
	               uvPrime2File << " " << uuPrime2LevelsList[hLevelsI].xy();
	               uwPrime2File << " " << uuPrime2LevelsList[hLevelsI].xz();
	               vvPrime2File << " " << uuPrime2LevelsList[hLevelsI].yy();
	               vwPrime2File << " " << uuPrime2LevelsList[hLevelsI].yz();
	               wwPrime2File << " " << uuPrime2LevelsList[hLevelsI].zz();

	               upPrime2File << " " << upPrime2LevelsList[hLevelsI].x();
	               vpPrime2File << " " << upPrime2LevelsList[hLevelsI].y();
	               wpPrime2File << " " << upPrime2LevelsList[hLevelsI].z();

	               uuuPrime3File << " " << uuuPrime3LevelsList[hLevelsI].xx();
	               uuvPrime3File << " " << uuuPrime3LevelsList[hLevelsI].xy();
	               uuwPrime3File << " " << uuuPrime3LevelsList[hLevelsI].xz();
	               uvvPrime3File << " " << uuuPrime3LevelsList[hLevelsI].yy();
	               uvwPrime3File << " " << uuuPrime3LevelsList[hLevelsI].yz();
	               uwwPrime3File << " " << uuuPrime3LevelsList[hLevelsI].zz();

	               vuuPrime3File << " " << vuuPrime3LevelsList[hLevelsI].xx();
	               vuvPrime3File << " " << vuuPrime3LevelsList[hLevelsI].xy();
	               vuwPrime3File << " " << vuuPrime3LevelsList[hLevelsI].xz();
	               vvvPrime3File << " " << vuuPrime3LevelsList[hLevelsI].yy();
	               vvwPrime3File << " " << vuuPrime3LevelsList[hLevelsI].yz();
	               vwwPrime3File << " " << vuuPrime3LevelsList[hLevelsI].zz();

	               wuuPrime3File << " " << wuuPrime3LevelsList[hLevelsI].xx();
	               wuvPrime3File << " " << wuuPrime3LevelsList[hLevelsI].xy();
	               wuwPrime3File << " " << wuuPrime3LevelsList[hLevelsI].xz();
	               wvvPrime3File << " " << wuuPrime3LevelsList[hLevelsI].yy();
	               wvwPrime3File << " " << wuuPrime3LevelsList[hLevelsI].yz();
	               wwwPrime3File << " " << wuuPrime3LevelsList[hLevelsI].zz();
	          }


                  UMeanFile << endl;
                  VMeanFile << endl;
                  WMeanFile << endl;
                  pMeanFile << endl;

                  uRMSFile << endl;
                  vRMSFile << endl;
                  wRMSFile << endl;
                  pRMSFile << endl;

                  uuPrime2File << endl;
                  uvPrime2File << endl;
                  uwPrime2File << endl;
                  vvPrime2File << endl;
                  vwPrime2File << endl;
                  wwPrime2File << endl;

                  upPrime2File << endl;
                  vpPrime2File << endl;
                  wpPrime2File << endl;

                  uuuPrime3File << endl;
                  uuvPrime3File << endl;
                  uuwPrime3File << endl;
                  uvvPrime3File << endl;
                  uvwPrime3File << endl;
                  uwwPrime3File << endl;

                  vuuPrime3File << endl;
                  vuvPrime3File << endl;
                  vuwPrime3File << endl;
                  vvvPrime3File << endl;
                  vvwPrime3File << endl;
                  vwwPrime3File << endl;

                  wuuPrime3File << endl;
                  wuvPrime3File << endl;
                  wuwPrime3File << endl;
                  wvvPrime3File << endl;
                  wvwPrime3File << endl;
                  wwwPrime3File << endl;
	     }
        }
