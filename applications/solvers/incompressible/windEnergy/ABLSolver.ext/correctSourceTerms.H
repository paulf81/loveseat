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
          
           SourceU = s;
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
               forAll(hLevelsCellList[level],i)
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
   }

   // This part is if the temperature source terms are applied directly as given
   if (temperatureSourceType == "given")
   {
       // This is if the source is only given at one height.  In that case,
       // assume the source is set uniformly throughout the domain to the
       // given value as a function of time only.
       if (nSourceTemperatureHeights == 1)
       {
           scalar sourceUX = interpolate2D(t,
                                           sourceHeightsTemperatureSpecified[0],
                                           sourceTemperatureTimesSpecified,
                                           sourceHeightsTemperatureSpecified,
                                           sourceTemperatureSpecified);
         
           SourceT = sourceT;
       }
       //  Otherwise, set the source as a function of height and time.
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
               forAll(hLevelsCellList[level],i)
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
