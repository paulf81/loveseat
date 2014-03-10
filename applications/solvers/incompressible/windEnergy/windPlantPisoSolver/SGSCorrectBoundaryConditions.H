//  Specify the SGS stress tensor and temperature flux vector at the lower surface using
//  Moeng's or Schummann & Grotzbach's model.
{
    const label patchID = mesh.boundaryMesh().findPatchID(lowerBoundaryName);
    const fvPatch& cPatch = mesh.boundary()[patchID];

    //  Get perpendicular distance from cell center to boundary.  In other words,
    //  the height of the z_1/2 grid level.
    dimensionedScalar z12
    (
         "z12", 
         dimensionSet(0,1,0,0,0,0,0), 
         hLevelsValues[0]
    );

    //  Get face normal vectors
    const vectorField normal = cPatch.nf();
    const vectorField normalArea = cPatch.Sf();

    //  Get face areas (individual and global sum) (note that "gSum" is used as opposed
    //  to "sum" because "gSum" is parallel-aware --- it gathers sums from each processor
    //  to which this patch belongs and makes the global sum)
    const scalarField area = cPatch.magSf();
    scalar areaTotal = gSum(area);

    //  Get resolved U vector at boundary face (UFace refers to patch().lookupPatchField...)
    const fvPatchVectorField& UFace =  cPatch.lookupPatchField<volVectorField, vector>("U");

    //  Get the components of the resolved velocity at z_1/2, U(z_1/2), that are
    //  parallel to the boundary face, U_||(z_1/2), by using:
    //
    //        U_||(z_1/2) = U(z_1/2) - (U(z_1/2) dot |S_n|)*|S_n|
    //
    //  where S_n is the surface normal vector.
    vectorField UParallel12 = UFace.patchInternalField();
    forAll(UFace, faceI)
    {
        UParallel12[faceI] = UParallel12[faceI] - ((UParallel12[faceI] & normal[faceI]) * normal[faceI]);
    }

    //  Get the magnitude of the U_||(z_1/2) field, its mean, and its fluctuation from the mean.
    scalarField magUParallel12 = mag(UParallel12);
    scalar meanMagUParallel12 = gSum(magUParallel12 * area)/areaTotal;
    scalarField fluctMagUParallel12 = magUParallel12;
    forAll(fluctMagUParallel12, faceI)
    {
        fluctMagUParallel12[faceI] = magUParallel12[faceI] - meanMagUParallel12;
    }
    scalar meanSqrFluctMagUParallel12 = gSum(Foam::sqr(fluctMagUParallel12) * area)/areaTotal;

    //  Find the mean over the surface of the U_||(z_1/2) vector and the
    //  magnitude of the mean vector (not the mean of the local magnitudes!!!)
    vector meanUParallel12 = gSum(UParallel12 * area)/areaTotal;
    dimensionedScalar magMeanUParallel12
    (
         "magMeanUParallel12",
         dimensionSet(0,1,-1,0,0,0,0),
         mag(meanUParallel12)
    );
    vectorField fluctUParallel12 = UParallel12;
    vectorField sqrFluctUParallel12 = UParallel12;
    forAll(fluctUParallel12, faceI)
    {
       fluctUParallel12[faceI] = UParallel12[faceI] - meanUParallel12;
       sqrFluctUParallel12[faceI].x() = fluctUParallel12[faceI].x()*fluctUParallel12[faceI].x();
       sqrFluctUParallel12[faceI].y() = fluctUParallel12[faceI].y()*fluctUParallel12[faceI].y();
       sqrFluctUParallel12[faceI].z() = fluctUParallel12[faceI].z()*fluctUParallel12[faceI].z();
    }
    vector meanSqrFluctUParallel12 = gSum(sqrFluctUParallel12 * area)/areaTotal;

    //  Calculate friction velocity, u*, using Monin-Obukhov similarity.  Note,
    //  that Monin-Obuhov assumes homogeneous z_1/2 and z_0.
    {
        const scalar eps = 0.1;
        const scalar tol = 1.0E-8;
        const int iterMax = 1000;
        int iter = 0;

        scalar L0 = 0.0;
        scalar L1 = 0.0;

        scalar phiM1_0 = 1.0;
        scalar phiM0_0 = 1.0;
        scalar psiM1_0 = 1.0;
        scalar psiM0_0 = 1.0;
        scalar phiM1_1 = 1.0;
        scalar phiM0_1 = 1.0;
        scalar psiM1_1 = 1.0;
        scalar psiM0_1 = 1.0;

        scalar f0 = 1.0E10 * tol;
        scalar f1 = 1.0E10 * tol;

        scalar gp0 = 0.0;

        scalar uStar0 = ((kappa * magMeanUParallel12) / Foam::log(z12 / z0)).value();
        scalar uStar1 = ((1.0 + eps) * uStar0);
        scalar uStar1Old = uStar1;

        dimensionedScalar uL
        (
             "uL",
             dimensionSet(0,1,0,0,0,0,0),
             1.0
        );
        dimensionedScalar uV
        (
             "uV",
             dimensionSet(0,1,-1,0,0,0,0),
             1.0
        );

        if (q0.value() == 0.0)       // neutral case
        {
            phiM1 = 1.0;
            uStar = (kappa * magMeanUParallel12) / Foam::log(z12 / z0);
        }

        else if (q0.value() > 0.0)   // unstable case
        {
            // Use a nonlinear Newton solver to solve for u*, L, phiM
            do
            {
                iter = iter + 1;

                L0      = -(Foam::pow(uStar0,3.0) / (kappa * (mag(g)/TRef) * q0)).value();
                phiM1_0 = Foam::pow((1.0 - (betaM * (z12 / L0).value())), -0.25);
                phiM0_0 = Foam::pow((1.0 - (betaM * (z0  / L0).value())), -0.25);
                psiM1_0 = 2.0 * Foam::log((1.0 + Foam::pow(phiM1_0, -1.0)) / 2.0) + Foam::log((1.0 + Foam::pow(phiM1_0, -2.0)) / 2.0) - 2.0 * Foam::atan(Foam::pow(phiM1_0, -1.0)) + Foam::constant::mathematical::pi/2.0;
                psiM0_0 = 2.0 * Foam::log((1.0 + Foam::pow(phiM0_0, -1.0)) / 2.0) + Foam::log((1.0 + Foam::pow(phiM0_0, -2.0)) / 2.0) - 2.0 * Foam::atan(Foam::pow(phiM0_0, -1.0)) + Foam::constant::mathematical::pi/2.0;
                f0      = magMeanUParallel12.value() - (uStar0/kappa) * ((log(z12 / z0)).value() - psiM1_0 + psiM0_0);

                L1      = -(Foam::pow(uStar1,3.0) / (kappa * (mag(g)/TRef) * q0)).value();
                phiM1_1 = Foam::pow((1.0 - (betaM * (z12 / L1).value())), -0.25);
                phiM0_1 = Foam::pow((1.0 - (betaM * (z0  / L1).value())), -0.25);
                psiM1_1 = 2.0 * Foam::log((1.0 + Foam::pow(phiM1_1, -1.0)) / 2.0) + Foam::log((1.0 + Foam::pow(phiM1_1, -2.0)) / 2.0) - 2.0 * Foam::atan(Foam::pow(phiM1_1, -1.0)) + Foam::constant::mathematical::pi/2.0;
                psiM0_1 = 2.0 * Foam::log((1.0 + Foam::pow(phiM0_1, -1.0)) / 2.0) + Foam::log((1.0 + Foam::pow(phiM0_1, -2.0)) / 2.0) - 2.0 * Foam::atan(Foam::pow(phiM0_1, -1.0)) + Foam::constant::mathematical::pi/2.0;
                f1      = magMeanUParallel12.value() - (uStar1/kappa) * ((log(z12 / z0)).value() - psiM1_1 + psiM0_1);

                gp0 = (f1 - f0) / (uStar1 - uStar0);

                uStar1Old = uStar1;
                uStar1 = uStar0 - (f0 / gp0);
                uStar0 = uStar1Old;

                uStar = uV * uStar1;
                L = uL * L1;
                phiM1 = phiM1_1;
            } while (mag(f1) > tol && iter < iterMax);

            if (iter == iterMax)
            {
                Info << "Max uStar iterations reached!!!" << endl;
            }
        }
        else if (q0.value() < 0.0)   // stable case
        {
            // Use a nonlinear Newton solver to solve for u*, L, phiM
            do
            {
                iter = iter + 1;

                L0      = -(Foam::pow(uStar0,3.0) / (kappa * (mag(g)/TRef) * q0)).value();
                phiM1_0 = 1.0 + (gammaM * (z12.value() / L0));
                f0      = magMeanUParallel12.value() - (uStar0/kappa) * ((log(z12 / z0)).value() + psiM1_0 - 1.0);

                L1      = -(Foam::pow(uStar1,3.0) / (kappa * (mag(g)/TRef) * q0)).value();
                phiM1_1 = 1.0 + (gammaM * (z12.value() / L1));
                f1      = magMeanUParallel12.value() - (uStar1/kappa) * ((log(z12 / z0)).value() + psiM1_1 - 1.0);

                gp0 = (f1 - f0) / (uStar1 - uStar0);

                uStar1Old = uStar1;
                uStar1 = uStar0 - (f0 / gp0);
                uStar0 = uStar1Old;

                uStar = uV * uStar1;
                L = uL * L1;
                phiM1 = phiM1_1;
            } while (mag(f1) > tol && iter < iterMax);

            if (iter == iterMax)
            {
                Info << "Max uStar iterations reached!!!" << endl;
            }
        }
    }


    //  Get resolved temperature at boundary face (TFace refers to patch().lookupPatchField...)
    const fvPatchScalarField& TFace =  cPatch.lookupPatchField<volScalarField, scalar>("T");

    //  Get the resolved temperature at z_1/2, T(z_1/2).
    scalarField T12 = TFace.patchInternalField();

    //  Find the mean over the surface of T(z_1/2)
    scalar meanT12 = gSum(T12 * area)/areaTotal;


    //  Apply the surface stress model.
    //  --Schummann & Grotzbach
    if(surfaceStressModel == "Schummann")
    {
        forAll(cPatch, faceI)
        {    
            devR.boundaryField()[patchID][faceI].xx() = 0.0;
            devR.boundaryField()[patchID][faceI].xy() = 0.0;
            devR.boundaryField()[patchID][faceI].xz() = -( (sqr(uStar)/magMeanUParallel12) * UParallel12[faceI].x() ).value();
            devR.boundaryField()[patchID][faceI].yy() = 0.0;
            devR.boundaryField()[patchID][faceI].yz() = -( (sqr(uStar)/magMeanUParallel12) * UParallel12[faceI].y() ).value();
            devR.boundaryField()[patchID][faceI].zz() = 0.0;
            q.boundaryField()[patchID][faceI].x()     = 0.0;
            q.boundaryField()[patchID][faceI].y()     = 0.0;
            q.boundaryField()[patchID][faceI].z()     = ( (q0/meanT12) * (T12[faceI]) ).value();
        }
    }
   
    //  --Moeng (1984)
    if(surfaceStressModel == "Moeng")
    {
        forAll(cPatch, faceI)
        {
            devR.boundaryField()[patchID][faceI].xx() = 0.0;
            devR.boundaryField()[patchID][faceI].xy() = 0.0;
            devR.boundaryField()[patchID][faceI].xz() = -( (sqr(uStar)/(meanMagUParallel12*magMeanUParallel12)) * (magUParallel12[faceI]*meanUParallel12.x()+meanMagUParallel12*(UParallel12[faceI].x()-meanUParallel12.x())) ).value();
            devR.boundaryField()[patchID][faceI].yy() = 0.0;
            devR.boundaryField()[patchID][faceI].yz() = -( (sqr(uStar)/(meanMagUParallel12*magMeanUParallel12)) * (magUParallel12[faceI]*meanUParallel12.y()+meanMagUParallel12*(UParallel12[faceI].y()-meanUParallel12.y())) ).value();
            devR.boundaryField()[patchID][faceI].zz() = 0.0;
            q.boundaryField()[patchID][faceI].x()     = 0.0;
            q.boundaryField()[patchID][faceI].y()     = 0.0;
            q.boundaryField()[patchID][faceI].z()     = ( (q0/(meanMagUParallel12*meanT12)) * (magUParallel12[faceI]*meanT12+meanMagUParallel12*(T12[faceI]-meanT12)) ).value();
        }
    }
 
    //  --Moeng (1984) with Wei/Brasseur (2010) correction
    if(surfaceStressModel == "MoengCorrected")
    {
        forAll(cPatch, faceI)
        {
            devR.boundaryField()[patchID][faceI].xx() = 0.0;
            devR.boundaryField()[patchID][faceI].xy() = 0.0;
            devR.boundaryField()[patchID][faceI].xz() = -( (sqr(uStar)/magMeanUParallel12) * meanUParallel12.x() ).value() 
                                                        +(  sqr(uStar) * betaSurfaceStress * ( -( fluctMagUParallel12[faceI]*meanUParallel12.x() + meanMagUParallel12*fluctUParallel12[faceI].x() )/
                                                                                                ( Foam::sqrt(meanSqrFluctMagUParallel12*Foam::sqr(meanUParallel12.x()) 
                                                                                                                           + Foam::sqr(meanMagUParallel12)*meanSqrFluctUParallel12.x()) ) ) ).value(); 
            devR.boundaryField()[patchID][faceI].yy() = 0.0;
            devR.boundaryField()[patchID][faceI].yz() = -( (sqr(uStar)/magMeanUParallel12) * meanUParallel12.y() ).value() 
                                                        +(  sqr(uStar) * betaSurfaceStress * ( -( fluctMagUParallel12[faceI]*meanUParallel12.y() + meanMagUParallel12*fluctUParallel12[faceI].y() )/
                                                                                                ( Foam::sqrt(meanSqrFluctMagUParallel12*Foam::sqr(meanUParallel12.y()) 
                                                                                                                           + Foam::sqr(meanMagUParallel12)*meanSqrFluctUParallel12.y()) ) ) ).value(); 
            devR.boundaryField()[patchID][faceI].zz() = 0.0;
            q.boundaryField()[patchID][faceI].x()     = 0.0;
            q.boundaryField()[patchID][faceI].y()     = 0.0;
            q.boundaryField()[patchID][faceI].z()     = ( (q0/(meanMagUParallel12*meanT12)) * (magUParallel12[faceI]*meanT12+meanMagUParallel12*(T12[faceI]-meanT12)) ).value();
        }
    }


    // Calculate some information to write to screen.
    symmTensor meanDevR = gSum(devR.boundaryField()[patchID] * area)/areaTotal;
    vector meanq = gSum(q.boundaryField()[patchID] * area)/areaTotal;
    Info << "C_s = " << Cs << endl;
    Info << "u* = " << uStar.value() << " m/s" << tab << "|<devR>|^(1/2) = " << Foam::pow(Foam::pow(meanDevR.xz(),2)+Foam::pow(meanDevR.yz(),2)+Foam::pow(meanDevR.zz(),2),0.25) << " m/s" << tab <<
            "<U>_1/2 = " << meanUParallel12 << " m/s" <<  endl;
    Info << "Qs = " << q0.value() << " K*(m/s)" << tab << "|<q>| = " << meanq.z() << " K*(m/s)" << tab <<
            "L = " << L.value() << " m"  << tab << "phi_m = " << phiM1  << tab << "<T>_1/2 = " << meanT12 << " K" << endl;
    
}




//  Specify the SGS stress tensor at the upper surface to be zero.
{
    const label patchID = mesh.boundaryMesh().findPatchID(upperBoundaryName);
    const fvPatch& cPatch = mesh.boundary()[patchID];
    forAll(cPatch, faceI)
    {
         devR.boundaryField()[patchID][faceI] = symmTensor::zero;	  

         q.boundaryField()[patchID][faceI] = vector::zero;
    }
}
