    // Declare the LES filter width variable.  Since stresses are computed on cell faces,
    // this needs to be on all cell faces.
    surfaceScalarField deltaLES
    (
        IOobject
        (
            "deltaLES",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("deltaLES", dimLength, SMALL)
    );

    {
        // Make a temporary LES filter width field at cell centers.  These are the filter
        // widths that are the cube-root of the volume of each cell.
        volScalarField deltaLESvol
        (
            IOobject
            (
                "deltaLESvol",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("deltaLESvol", dimLength, SMALL)
        );

        // Do the cube root.
        deltaLESvol.internalField() = pow(mesh.V(),1.0/3.0);

        // Update boundary values (mainly to make sure processor boundaries are updated).
        deltaLESvol.correctBoundaryConditions();

        // Linearly interpolate the cell center LES filter widths to the cell faces.
        deltaLES = linearInterpolate(deltaLESvol);

        // Now take care of boundary values.  Cyclic and processor boundaries will be
        // interpolated properly.  Other types of boundaries don't know what to do, so
        // they default to a value of "SMALL" = 1E-15.  Take the value from the owner
        // cell as the value at the face.
        forAll(mesh.boundary(),patchID)
        {
            const fvPatch& cPatch = mesh.boundary()[patchID];
            const unallocLabelList& faceCellsPatch = cPatch.faceCells();

	    if (cPatch.coupled() == false)
            {
                forAll(cPatch,faceI)
                {
                    label cellI = faceCellsPatch[faceI];
                    deltaLES.boundaryField()[patchID][faceI] = deltaLESvol[cellI];
                }
            }
        }

        // Multiply by a filter width scaling factor.
        deltaLES *= deltaLESCoeff;
    }
