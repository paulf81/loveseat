    // This piece of code finds the vertical faces that are adjacent to the
    // upper and lower boundaries. It is for use later in correcting the SGS
    // viscosity and temperature diffusivity at this level since direct
    // calculation of velocity and temperature gradients (required in 
    // computing SGS viscosity and temperature diffusivity) is not
    // accurately possible (i.e. u,v, and T are not required at surface
    // in order to compute momentum and temperature transport equations, so
    // there is no way to get gradient at first level (z = 1/2) away from
    // surface.)

    // The adjacent vertical faces are divided into two groups -- those
    // adjacent to the lower and upper boundaries.  Within each group are
    // two subgroups -- those faces that are internal to the mesh and those
    // that lie on a mesh boundary.

    // Bottom adjacent faces.
    //   internal faces.
    DynamicList<label> bottomAdjInternalFaces(0);
    forAll(mesh.owner(), faceI)
    {
         if(mag(mesh.Cf()[faceI][upIndex] - hLevelsValues[0]) <= hLevelsTol)
         {
              bottomAdjInternalFaces.append(faceI);
         }

    }
  
    //   boundary faces.
    List<DynamicList<label> > bottomAdjBoundaryFaces(mesh.boundaryMesh().size());
    List<label> bottomAdjBoundaryFacesCount(mesh.boundaryMesh().size());
    forAll(mesh.boundaryMesh(), patchID)
    {
         const fvPatch& cPatch = mesh.boundary()[patchID];
         bottomAdjBoundaryFacesCount[patchID] = 0;
         forAll(cPatch, faceI)
         {
              if(mag(cPatch.Cf()[faceI][upIndex] - hLevelsValues[0]) <= hLevelsTol)
              {
                   bottomAdjBoundaryFaces[patchID].append(faceI);
                   bottomAdjBoundaryFacesCount[patchID]++;
              }
         }
    }


    // Top adjacent faces.
    //   internal faces.
    DynamicList<label> topAdjInternalFaces(0) ;
    forAll(mesh.owner(), faceI)
    {
         if(mag(mesh.Cf()[faceI][upIndex] - hLevelsValues[hLevelsTotal-1]) <= hLevelsTol)
         {
              topAdjInternalFaces.append(faceI) ;
         }

    }

    //   boundary faces.
    List<DynamicList<label> > topAdjBoundaryFaces(mesh.boundaryMesh().size());
    List<label> topAdjBoundaryFacesCount(mesh.boundaryMesh().size());
    forAll(mesh.boundaryMesh(), patchID)
    {
         const fvPatch& cPatch = mesh.boundary()[patchID];
         topAdjBoundaryFacesCount[patchID] = 0;
         forAll(cPatch, faceI)
         {
              if(mag(cPatch.Cf()[faceI][upIndex] - hLevelsValues[hLevelsTotal - 1]) <= hLevelsTol)
              {
                   topAdjBoundaryFaces[patchID].append(faceI);
                   topAdjBoundaryFacesCount[patchID]++;
              }
         }
    }


    // Write information to log file or screen for user verification.
    //Pout << "Total Mesh Faces Count = " << mesh.faces().size() << endl << endl;
    //Pout << "Bottom Adjacent Internal Face Count = " << bottomAdjInternalFaces.size() << endl;
    label bottomAdjTotalFacesCount = bottomAdjInternalFaces.size();
    forAll(mesh.boundaryMesh(), patchID)
    {
         //Pout << "Bottom Adjacent Boundary [" << patchID << "] Face Count = " << bottomAdjBoundaryFacesCount[patchID] << endl;
         bottomAdjTotalFacesCount += bottomAdjBoundaryFacesCount[patchID];
    }
    //Pout << "Bottom Adjacent Total Face Count = " << bottomAdjTotalFacesCount << endl << endl;
    //Pout << "Top Adjacent Internal Face Count    = " << topAdjInternalFaces.size() << endl;
    label topAdjTotalFacesCount = topAdjInternalFaces.size();
    forAll(mesh.boundaryMesh(), patchID)
    {
         //Pout << "Top Adjacent Boundary [" << patchID << "] Face Count    = " << topAdjBoundaryFacesCount[patchID] << endl;
         topAdjTotalFacesCount += topAdjBoundaryFacesCount[patchID];
    }
    //Pout << "Top Adjacent Total Face Count    = " << topAdjTotalFacesCount << endl << endl;

