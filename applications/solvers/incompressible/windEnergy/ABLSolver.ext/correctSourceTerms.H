{
   // Get the current time.
   scalar t = runTime.value();

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
       }
       //  Otherwise, set the source as a function of height and time.
       else
       {
          // Interpolate the source values in time and height.
           List<scalar> sourceUXColumn = interpolate2D(t,
                                                       hLevelsValues,
                                                       sourceMomentumXTimesSpecified,
                                                       sourceHeightsMomentumSpecified,
                                                       sourceMomentumXSpecified);
 
           List<scalar> sourceUYColumn = interpolate2D(t,
                                                       hLevelsValues,
                                                       sourceMomentumYTimesSpecified,
                                                       sourceHeightsMomentumSpecified,
                                                       sourceMomentumYSpecified);

           List<scalar> sourceUZColumn = interpolate2D(t, 
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

               // add the correction on to the source field.
               for (label i = 0; i < numCellPerLevel[level]; i++)
               {
                   dsVol[hLevelsCellList[level][i]] = ds;
                   SourceU[hLevelsCellList[level][i]] += ds;
               }
           }
           dsVol.correctBoundaryConditions();
           SourceU.correctBoundaryConditions();
          // if (Pstream::myProcNo() == 0)
          // {
          //    Pout << dsVol << endl;
          // }
           

           // Explicitly correct velocity and flux.
           U += rAU * dsVol;
           U.correctBoundaryConditions();
           

//    surfaceScalarField dsVolFlux = rAUf * (fvc::interpolate(dsVol) & mesh.Sf());
//    dsVolFlux.correctBoundaryConditions();
//    scalar aa;
//    forAll(dsVolFlux.boundaryField(), patchi)
//    {
//        scalar sumPhiBoundary = 0.0;
//        const fvsPatchScalarField& phip = dsVolFlux.boundaryField()[patchi];
//        forAll(phip,i)
//        {
//           sumPhiBoundary += phip[i];
//        }
//        aa += sumPhiBoundary;
//    }

  //  Pout << "Boundary dsVol Flux on Processor " << Pstream::myProcNo() << " = " << aa << endl;

          // if (Pstream::myProcNo() == 0)
          // {
          //    Pout << "interp = " << rAUf * (fvc::interpolate(dsVol) & mesh.Sf());
          // }
           phi += rAUf * (fvc::interpolate(dsVol) & mesh.Sf());
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
       }
       // Otherwise, set the source as a function of height and time.
       else
       {
          // Interpolate the source values in time and height.
           List<scalar> sourceTColumn = interpolate2D(t,
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
       }
   }

   // This part is if the temperature source terms have to be computed.
   else if (temperatureSourceType == "computed")
   {
  
   }
}
