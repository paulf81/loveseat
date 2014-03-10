// Calculate divergence of velocity and velocity flux and display
{
    divPhi = fvc::div(phi);
    volScalarField divPhiMag  = pow(pow(divPhi,2),0.5);
    scalar minLocalPhiContErr = min(divPhiMag).value();
    reduce(minLocalPhiContErr, minOp<scalar>());
    scalar maxLocalPhiContErr = max(divPhiMag).value();
    reduce(maxLocalPhiContErr, maxOp<scalar>());
    scalar avgLocalPhiContErr = divPhiMag.weightedAverage(mesh.V()).value();
    Info << "Local Flux Continuity Error:  Min " << minLocalPhiContErr << tab
         <<                               "Max " << maxLocalPhiContErr << tab
         <<                     "Weighted Mean " << avgLocalPhiContErr << endl;

    divU   = fvc::div(U);
    volScalarField divUMag    = pow(pow(divU,2),0.5);
    scalar minLocalUContErr   = min(divUMag).value();
    reduce(minLocalUContErr, minOp<scalar>());
    scalar maxLocalUContErr   = max(divUMag).value();
    reduce(maxLocalUContErr, maxOp<scalar>());
    scalar avgLocalUContErr   = divUMag.weightedAverage(mesh.V()).value();
    Info << "Local U Continuity Error:     Min " << minLocalUContErr   << tab
         <<                               "Max " << maxLocalUContErr << tab
         <<                     "Weighted Mean " << avgLocalUContErr   << endl;

    scalar globalSumPhiBoundary = 0.0;
    forAll(phi.boundaryField(), patchi)
    {
        scalar sumPhiBoundary = 0.0;
        const fvsPatchScalarField& phip = phi.boundaryField()[patchi];
        forAll(phip,i)
        {
           sumPhiBoundary += phip[i];
        }
        globalSumPhiBoundary += sumPhiBoundary;
      //Pout << "Boundary " << mesh.boundaryMesh()[patchi].name() << " Flux:  " << sumPhiBoundary << endl;
    }

    reduce(globalSumPhiBoundary, sumOp<scalar>());
    Info << "Total Boundary Flux: " << globalSumPhiBoundary << endl;
}
