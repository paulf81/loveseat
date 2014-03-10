//  Update the SGS momentum and temperature fluxes using the Smagorinsky model.
{
    // Calculate the velocity gradient tensor on cell faces.
    surfaceTensorField gradU = (surfNorm * fvc::snGrad(U))  + ((identity - surfNorm*surfNorm) & (fvc::interpolate(fvc::grad(U))));

    // Calculate the temperature gradient on cell faces.
    surfaceVectorField gradT = (surfNorm * fvc::snGrad(T))  + ((identity - surfNorm*surfNorm) & (fvc::interpolate(fvc::grad(T))));



    // The SGS viscosity at the center of the vertical faces adjacent to the bottom boundary
    // (grid level z_1/2) will not be correct since the velocity gradient calculated there 
    // is fairly meaningless.  It is meaningless since there is no need to specify horizontal 
    // velocity at the surface (it is not no slip since it is a rough wall), so one cannot 
    // get an accurate estimate of the velocity gradient at this level.  
    //
    // One way to deal with this is to set the velocity gradient at the z_1/2 level equal to
    // that at the z_1 (top of first cell adjacent to surface) level, meaning that the
    // SGS viscosity would be the same at z_1/2 and z_1.  Or, one could just set the 
    // SGS viscosity z_1/2 to be the same as that at z_1, which is what is done here.  The
    // SGS viscosity points at z_1 are different horizontally than those at z_1/2, so
    // linear interpolation is necessary.
    //
    // One last idea, which is not implemented here, but which may be a good one would be
    // to use Monin-Obuhkov similarity locally to calculate the velocity gradient at
    // the z_1/2 level (similar to Moeng, 1984, pg. 2056, section 5, first paragraph).
    //
    // The same is done for the SGS temperature diffusivity.
    //
    // Here is the process to do this:
    {
        // Find the patch that corresponds to the lower boundary.
        label patchID = mesh.boundaryMesh().findPatchID(lowerBoundaryName);
        const fvPatch& cPatch = mesh.boundary()[patchID];

        // Copy the vertical components of the velocity gradient tensor and the temperature
        // gradient vector from the z_1 level to the z_1/2 cell-center level.  Temporarily 
        // store these in the divDevR vector and the divq scalar, respectively, at the z_1/2
        // cell-center level.
        forAll(cPatch, faceI)
        {
            // Get the label for the cell to which this patch face belongs.
            label cellI = cPatch.faceCells()[faceI];

            // Get the label for the face opposite this patch face.
            label oppFaceI = mesh.cells()[cellI].opposingFaceLabel(faceI+cPatch.patch().start(),mesh.faces());

            // Get the part of the velocity gradient that is on the opposite face (at the z_1 level) that is
            // normal to this patch face.
            // Sometimes the opposite face is a patch face (in places where the subdomain is only one
            // layer of cells thick.  In that case, the opposingFaceLabel returns the absolute face label
            // for the opposite face, but we need the label local to the patch that it lies on and the label
            // of that patch.  Then we use boundaryField() to access the gradient there.  Otherwise, we don't
            // need to do all of this.
            if(oppFaceI >= gradU.size())
            {
                // Get the label for the patch that the opposite face lies upon.
                label oppPatchID = mesh.boundaryMesh().whichPatch(oppFaceI);

                // Get the local label for the opposite face.
                oppFaceI -= mesh.boundary()[oppPatchID].patch().start();

                // Then copy the z_1 level vertical velocity gradient into divDevR.
                divDevR[cellI] = -(surfNorm.boundaryField()[patchID][faceI] & gradU.boundaryField()[oppPatchID][oppFaceI]);
                divq[cellI]    = -(surfNorm.boundaryField()[patchID][faceI] & gradT.boundaryField()[oppPatchID][oppFaceI]);
            }
            else
            {
                divDevR[cellI] = -(surfNorm.boundaryField()[patchID][faceI] & gradU[oppFaceI]);
                divq[cellI] = -(surfNorm.boundaryField()[patchID][faceI] & gradT[oppFaceI]);
            }

        }
        // We have modified the internal field of divDevR and divq.  If we call interpolate
        // on these field, OpenFOAM will interpolate to interior boundaries just fine, but it 
        // will not interpolate to outer boundaries (like processor boundaries).  It will not
        // do this because the values there are from boundary conditions.  Those might be
        // "calculated" type boundary conditions that are interpolated from the interior, but
        // unless we explicitly update those boundary conditions, the boundary values (and
        // hence the interpolate command) will not take into account the update we made by
        // copying from the z_1 to z_1/2 level.
        divDevR.correctBoundaryConditions();
        divq.correctBoundaryConditions();

        // We need to interpolate the vertical components of the velocity gradient tensor 
        // that were copied to the z_1/2 cell-center level to the z_1/2 vertical face centers.
        // Do this by interpolating from cell centers to cell faces into temporary variables
        // defined at cell surfaces. (There is probably a more efficient way to do this instead
        // of interpolating everywhere!)
        surfaceVectorField gradUTmp = fvc::interpolate(divDevR);
        surfaceScalarField gradTTmp = fvc::interpolate(divq);

        // Then for all internal vertical wall-adjacent faces, copy from the temporary variable
        // into the real variables gradU and gradT.
        forAll(bottomAdjInternalFaces, j)
        {
            label faceI = bottomAdjInternalFaces[j];

            gradU[faceI].zx() = gradUTmp[faceI].x();
            gradU[faceI].zy() = gradUTmp[faceI].y();
            gradU[faceI].zz() = gradUTmp[faceI].z();

            gradT[faceI].z()  = gradTTmp[faceI];
        }  

        // Do the same at the boundary vertical wall-adjacent faces.
        forAll(mesh.boundaryMesh(), patchID)
        {
            forAll(bottomAdjBoundaryFaces[patchID], j)
            {
                label faceI = bottomAdjBoundaryFaces[patchID][j];

                gradU.boundaryField()[patchID][faceI].zx() = gradUTmp.boundaryField()[patchID][faceI].x();
                gradU.boundaryField()[patchID][faceI].zy() = gradUTmp.boundaryField()[patchID][faceI].y();
                gradU.boundaryField()[patchID][faceI].zz() = gradUTmp.boundaryField()[patchID][faceI].z();

                gradT.boundaryField()[patchID][faceI].z()  = gradTTmp.boundaryField()[patchID][faceI];
            }
        }
    }




    // Now that the gradients at the z_1 level are good, calculate the SGS viscosity throughout 
    // the domain on cell faces.
    nuLES = sqr(Cs*deltaLES)*sqrt(2.0*magSqr(dev(symm(gradU))));

    // Next, compute the turbulent Prandtl number throughout the domain on all cell faces.
    // This is a bit more involved that nuLES since it is sensitized to stability.
    {
         // first, compute turbulent Prandtl number (sensitized to stratification)
         // requiring an SGS velocity scale.  From Moeng (1984) one can equate the
         // velocity scale 0.1*sqrt(e') to the Smagorinsky velocity scale Cs^2*delta*||S_ij||.
         // Then that equality is used in Moeng's equation 18b to get the length
         // scale ls.
         surfaceScalarField velScale = nuLES/deltaLES;

         // also, a measure of stability is needed.
         dimensionedScalar tiny
         (
             "tiny",
             dimensionSet(0,0,-2,0,0,0,0),
             scalar(1.0E-20)
         );
         surfaceScalarField stability = (mag(g)/TRef) * (gradT & nUp);

         // and, a lengthscale is needed
         surfaceScalarField lPos = max(min(7.6*velScale*Foam::sqrt(1.0/max(tiny,stability)),deltaLES),0.0075*deltaLES);
         surfaceScalarField lNeg = 0.0*lPos + deltaLES;
         surfaceScalarField l = (1.0-pos(stability))*lNeg + pos(stability)*lPos;
         surfaceScalarField Prt = 1.0/(1.0 + 2.0*(l/deltaLES));
         
         // Then use the SGS viscosity and turbulent Prandtl number to compute the turbulent
         // diffusivity of temperature at cell faces.
         kappaLES = nuLES/Prt;
    }

    // Update the deviatoric part of the SGS stress tensor.
    surfaceSymmTensorField devR = -nuLES*dev(twoSymm(gradU));

    // Update the SGS temperature flux vector.
    surfaceVectorField q = -kappaLES*gradT;
       
    // Update the deviatoric part of the SGS stress tensor and the SGS temperature
    // flux vector at lower and upper faces.  The upper boundary stress and temperature
    // flux is set to zero.  The lower boundary uses a rough wall model.
    #include "SGSCorrectBoundaryConditions.H"

    // Calculate the divergence of the SGS stress tensor
    divDevR = fvc::div(mesh.Sf() & devR);

    // Calculate the divergence of the SGS temperature flux vector.
    divq = fvc::div(mesh.Sf() & q);



    // Just a check to see if fluxes are being assembled correctly.
    //{
    //    const unallocLabelList& owner = mesh.owner();
    //    const unallocLabelList& neighbour = mesh.neighbour();
    //    label cellLabel = 0;
    //    vector sum = vector::zero;
    //    surfaceVectorField devRDotSf = devR & mesh.Sf();
    //    forAll(devRDotSf,facei)
    //    {
    //        if (owner[facei] == cellLabel)
    //        {
    //            Info << "Interior owner " << facei << " " << devRDotSf[facei] << endl;
    //            sum += devRDotSf[facei];
    //        }
    //        if (neighbour[facei] == cellLabel)
    //        {
    //            Info << "Interior neighbour " << facei << " " << devRDotSf[facei] << endl;
    //            sum -= devRDotSf[facei];
    //        }
    //    }
    //    forAll(mesh.boundary(),patchi)
    //    {
    //        const unallocLabelList& pFaceCells = mesh.boundary()[patchi].faceCells();
    //        forAll(pFaceCells,facei)
    //        {
    //            if(pFaceCells[facei] == cellLabel)
    //            {
    //                Info << "Boundary " << patchi << " face " << facei << " " << devRDotSf.boundaryField()[patchi][facei] << endl;
    //                sum += devRDotSf.boundaryField()[patchi][facei];
    //            }
    //        }
    //    }
    //    Info << divDevR[cellLabel] << " " << sum/mesh.V()[cellLabel] << " " << divDevR[cellLabel] - sum/mesh.V()[cellLabel] << endl;
    //}
}	
