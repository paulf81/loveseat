{
   // Get the current time.
   scalar t = runTime.value();

   // This part is if the momentum source terms are applied directly as given
   if (momentumSourceType == "given")
   {
      if (nSourceMomentumHeights == 1)
      {
      }
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
   }

   // This part is if the temperature source terms have to be computed.
   else if (temperatureSourceType == "computed")
   {
   }
}
