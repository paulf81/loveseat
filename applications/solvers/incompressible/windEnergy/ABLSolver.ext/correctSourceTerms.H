{
   // Get the current time and time step size.
   scalar t = runTime.value();
   scalar dt = runTime.deltaT().value();

   // This part is if the momentum source terms are applied directly as given
   if (momentumSourceType == "given")
   {
       // This is if the source is only given at one height.  In that case,
       // assume the source is set uniformly throughout the domain to the
       // given value as a function of time only.
       if (nSourceMomentumHeights == 1)
       {
           scalar sourceUX = interpolate2D(t,
                                           sourceHeightsMomentumSpecified[0],
                                           sourceMomentumXTimesSpecified,
                                           sourceHeightsMomentumSpecified,
                                           sourceMomentumXSpecified);

           scalar sourceUY = interpolate2D(t,
                                           sourceHeightsMomentumSpecified[0],
                                           sourceMomentumYTimesSpecified,
                                           sourceHeightsMomentumSpecified,
                                           sourceMomentumYSpecified);

           scalar sourceUZ = interpolate2D(t,
                                           sourceHeightsMomentumSpecified[0],
                                           sourceMomentumZTimesSpecified,
                                           sourceHeightsMomentumSpecified,
                                           sourceMomentumZSpecified);
         
           vector s(vector::zero);
           s.x() = sourceUX;
           s.y() = sourceUY;
           s.z() = sourceUZ;

           forAll(SourceU,i)
           {
               SourceU[i] = s;
           }


           // Write the source information.
           if (Pstream::master())
           {
               if (statisticsOn)
               {
                   if (runTime.timeIndex() % statisticsFreq == 0)
                   {
                       sourceUXHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value() << " " << sourceUX << endl;
                       sourceUYHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value() << " " << sourceUY << endl;
                       sourceUZHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value() << " " << sourceUZ << endl;
                   }
               }
           }
       }
       //  Otherwise, set the source as a function of height and time.
       else
       {
          // Interpolate the source values in time and height.
           sourceUXColumn = interpolate2D(t,
                                          hLevelsValues,
                                          sourceMomentumXTimesSpecified,
                                          sourceHeightsMomentumSpecified,
                                          sourceMomentumXSpecified);
 
           sourceUYColumn = interpolate2D(t,
                                          hLevelsValues,
                                          sourceMomentumYTimesSpecified,
                                          sourceHeightsMomentumSpecified,
                                          sourceMomentumYSpecified);

           sourceUZColumn = interpolate2D(t, 
                                          hLevelsValues,
                                          sourceMomentumZTimesSpecified,
                                          sourceHeightsMomentumSpecified,
                                          sourceMomentumZSpecified);

           // Now go by cell levels and apply the source term.
           forAll(hLevelsCellList,level)
           {
               for (label i = 0; i < numCellPerLevel[level]; i++)
               {
                   vector s(vector::zero);
                   s.x() = sourceUXColumn[level];
                   s.y() = sourceUYColumn[level];
                   s.z() = sourceUZColumn[level];

                   SourceU[hLevelsCellList[level][i]] = s;
               }
           }                               


           // Write the column of source information.
           if (Pstream::master())
           {
               if (statisticsOn)
               {
                   if (runTime.timeIndex() % statisticsFreq == 0)
                   {
                       sourceUXHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value();
                       sourceUYHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value();
                       sourceUZHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value();

                       forAll(hLevelsCellList,level)
                       {
                          sourceUXHistoryFile() << " " << sourceUXColumn[level];
                          sourceUYHistoryFile() << " " << sourceUYColumn[level];
                          sourceUZHistoryFile() << " " << sourceUZColumn[level];
                       }

                       sourceUXHistoryFile() << endl;
                       sourceUYHistoryFile() << endl;
                       sourceUZHistoryFile() << endl;
                   }
               }
           }
       }
   }
   
   // This part is if the momentum source terms have to be computed.
   else if (momentumSourceType == "computed")
   {
       // This is if the velocity is only specified at one height.  In that case,
       // a source term will be computed and applied uniformly in all three
       // dimensions, the same as in the original ABLSolver.
       if (nSourceMomentumHeights == 1)
       {
       }
       // Otherwise, set the source as a function of height and time.
       else
       {
           // Get the inverse of the A operator on the U equation matrix.
           // The A operator gives the diagonal elements of the matrix.
           // Also, interpolate it to the cell faces.
           volScalarField rAU("rAU", 1.0/UEqn.A());
           surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

       
           // Create a volVectorField for the source term change.
           volVectorField dsVol("dsVol", SourceU);


           // Interpolate the desired velocity values in time and height.
           List<scalar> desiredUXColumn = interpolate2D(t,
                                                        hLevelsValues,
                                                        sourceMomentumXTimesSpecified,
                                                        sourceHeightsMomentumSpecified,
                                                        sourceMomentumXSpecified);

           List<scalar> desiredUYColumn = interpolate2D(t,
                                                        hLevelsValues,
                                                        sourceMomentumYTimesSpecified,
                                                        sourceHeightsMomentumSpecified,
                                                        sourceMomentumYSpecified);

           List<scalar> desiredUZColumn = interpolate2D(t,
                                                        hLevelsValues,
                                                        sourceMomentumZTimesSpecified,
                                                        sourceHeightsMomentumSpecified,
                                                        sourceMomentumZSpecified);

            
           // Compute the planar-averaged actual velocity and A operator at each cell level.
           List<vector> Umean(hLevelsTotal,vector::zero);
           List<vector> UmeanDesired(hLevelsTotal,vector::zero);
           List<scalar> rAUmean(hLevelsTotal,0.0);
   

           forAll(hLevelsCellList,level)
           {
               for (label i = 0; i < numCellPerLevel[level]; i++)
               {
                   Umean[level] += U[hLevelsCellList[level][i]] * mesh.V()[hLevelsCellList[level][i]];
                   rAUmean[level] += rAU[hLevelsCellList[level][i]] * mesh.V()[hLevelsCellList[level][i]];
               }
           }

           reduce(Umean,sumOp<List<vector> >());
           reduce(rAUmean,sumOp<List<scalar> >());

           forAll(hLevelsCellList,level)
           {
               Umean[level] /= totVolPerLevel[level];
               rAUmean[level] /= totVolPerLevel[level];

               vector v(vector::zero);
               v.x() = desiredUXColumn[level];
               v.y() = desiredUYColumn[level];
               v.z() = desiredUZColumn[level];
               UmeanDesired[level] = v;
           }

           Info << "UmeanDesired = " << UmeanDesired << endl;
           Info << "UmeanActual = " << Umean << endl;


           // Compute the correction to the source term.
           forAll(hLevelsCellList,level)
           {
               // this is the correction a this level.
               vector ds = (UmeanDesired[level] - Umean[level]) / rAUmean[level];

               // subtract off any vertical part.
               ds -= (ds & nUp) * nUp;

               // apply relaxation.
               ds *= alphaMomentum;

               // add the correction to the source column vector.
               sourceUXColumn[level] += ds.x();
               sourceUYColumn[level] += ds.y();
               sourceUZColumn[level] += ds.z();

               // add the correction on to the source field.
               for (label i = 0; i < numCellPerLevel[level]; i++)
               {
                   dsVol[hLevelsCellList[level][i]] = ds;
                   SourceU[hLevelsCellList[level][i]] += ds;
               }
           }
           dsVol.correctBoundaryConditions();
           SourceU.correctBoundaryConditions();

           // Explicitly correct velocity and flux.
           U += rAU * dsVol;
           U.correctBoundaryConditions();
           phi += rAUf * (fvc::interpolate(dsVol) & mesh.Sf());

           // Write the column of source information.
           if (Pstream::master())
           {
               if (statisticsOn)
               {
                   if (runTime.timeIndex() % statisticsFreq == 0)
                   {
                       sourceUXHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value();
                       sourceUYHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value();
                       sourceUZHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value();

                       forAll(hLevelsCellList,level)
                       {
                          sourceUXHistoryFile() << " " << sourceUXColumn[level];
                          sourceUYHistoryFile() << " " << sourceUYColumn[level];
                          sourceUZHistoryFile() << " " << sourceUZColumn[level];
                       }

                       sourceUXHistoryFile() << endl;
                       sourceUYHistoryFile() << endl;
                       sourceUZHistoryFile() << endl;
                   }
               }
           }
       }
   }


   // This part is if the temperature source terms are applied directly as given
   if (temperatureSourceType == "given")
   {
       // This is if the source is only given at one height.  In that case,
       // assume the source is set uniformly throughout the domain to the
       // given value as a function of time only.
       if (nSourceTemperatureHeights == 1)
       {
           scalar sourceT = interpolate2D(t,
                                          sourceHeightsTemperatureSpecified[0],
                                          sourceTemperatureTimesSpecified,
                                          sourceHeightsTemperatureSpecified,
                                          sourceTemperatureSpecified);
         
           forAll(SourceT,i)
           {           
               SourceT[i] = sourceT;
           }

           // Write the source information.
           if (Pstream::master())
           {
               if (statisticsOn)
               {
                   if (runTime.timeIndex() % statisticsFreq == 0)
                   {
                       sourceTHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value() << " " << sourceT << endl;
                   }
               }
           }
       }
       // Otherwise, set the source as a function of height and time.
       else
       {
           // Interpolate the source values in time and height.
           sourceTColumn = interpolate2D(t,
                                         hLevelsValues,
                                         sourceTemperatureTimesSpecified,
                                         sourceHeightsTemperatureSpecified,
                                         sourceTemperatureSpecified);

           // Now go by cell levels and apply the source term.
           forAll(hLevelsCellList,level)
           {
               for (label i = 0; i < numCellPerLevel[level]; i++)
               {
                   SourceT[hLevelsCellList[level][i]] = sourceTColumn[level];
               }
           }                               

           if (Pstream::master())
           {
               if (statisticsOn)
               {
                   if (runTime.timeIndex() % statisticsFreq == 0)
                   {
                       sourceTHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value();

                       forAll(hLevelsCellList,level)
                       {
                          sourceTHistoryFile() << " " << sourceTColumn[level];
                       }

                       sourceTHistoryFile() << endl;
                   }
               }
           }
       }
   }

   // This part is if the temperature source terms have to be computed.
   else if (temperatureSourceType == "computed")
   {
       // This is if the temperature is only specified at one height.  In that case,
       // a source term will be computed and applied uniformly in all three
       // dimensions.
       if (nSourceMomentumHeights == 1)
       {
       }
       // Otherwise, set the source as a function of height and time.
       else
       {
           // Create a volScalarField for the source term change.
           volScalarField dsVol("dsVol", SourceT);


           // Interpolate the desired velocity values in time and height.
           List<scalar> desiredTColumn = interpolate2D(t,
                                                       hLevelsValues,
                                                       sourceTemperatureTimesSpecified,
                                                       sourceHeightsTemperatureSpecified,
                                                       sourceTemperatureSpecified);

            
           // Compute the planar-averaged actual temperature at each cell level.
           List<scalar> Tmean(hLevelsTotal,0.0);
           List<scalar> TmeanDesired(hLevelsTotal,0.0);
   

           forAll(hLevelsCellList,level)
           {
               for (label i = 0; i < numCellPerLevel[level]; i++)
               {
                   Tmean[level] += T[hLevelsCellList[level][i]] * mesh.V()[hLevelsCellList[level][i]];
               }
           }

           reduce(Tmean,sumOp<List<scalar> >());

           forAll(hLevelsCellList,level)
           {
               Tmean[level] /= totVolPerLevel[level];

               TmeanDesired[level] = desiredTColumn[level];
           }

           Info << "TmeanDesired = " << TmeanDesired << endl;
           Info << "TmeanActual = " << Tmean << endl;

        
           // Compute the correction to the source term.
           forAll(hLevelsCellList,level)
           {
               // this is the correction a this level.
               scalar ds = (TmeanDesired[level] - Tmean[level]) / dt;

               // apply relaxation.
               ds *= alphaTemperature;

               // add the correction to the source column vector.
               sourceTColumn[level] += ds;

               // add the correction on to the source field.
               for (label i = 0; i < numCellPerLevel[level]; i++)
               {
                   dsVol[hLevelsCellList[level][i]] = ds;
                   SourceT[hLevelsCellList[level][i]] += ds;
               }
           }
           dsVol.correctBoundaryConditions();
           SourceT.correctBoundaryConditions();

           // Explicitly correct velocity and flux.
           T += runTime.deltaT() * dsVol;
           T.correctBoundaryConditions();

           // Write the column of source information.
           if (Pstream::master())
           {
               if (statisticsOn)
               {
                   if (runTime.timeIndex() % statisticsFreq == 0)
                   {
                       sourceTHistoryFile() << runTime.timeName() << " " << runTime.deltaT().value();

                       forAll(hLevelsCellList,level)
                       {
                          sourceTHistoryFile() << " " << sourceTColumn[level];
                       }

                       sourceTHistoryFile() << endl;
                   }
               }
           }
       }
   }
}
