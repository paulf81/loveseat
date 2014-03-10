/*---------------------------------------------------------------------------*\
This file was modified or created at the National Renewable Energy
Laboratory (NREL) on January 6, 2012 in creating the SOWFA (Simulator for
Offshore Wind Farm Applications) package of wind plant modeling tools that
are based on the OpenFOAM software. Access to and use of SOWFA imposes
obligations on the user, as set forth in the NWTC Design Codes DATA USE
DISCLAIMER AGREEMENT that can be found at
<http://wind.nrel.gov/designcodes/disclaimer.html>.
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    windPlantPisoSolver 

Description
    Transient solver for incompressible wind plant flows.

    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
//#include "turbulenceModel.H"
#include "LESModel.H"
#include "IFstream.H"
#include "OFstream.H"
#include "horizontalAxisWindTurbinesALM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalCoriolisAcceleration.H"
    #include "createFields.H"
    #include "createGradPd.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    // Get/compute initial information
    //   The identity matrix of rank 3.
//  dimensionedTensor identity
//  (
//       "identity", 
//       dimensionSet(0,0,0,0,0,0,0), 
//       tensor(1,0,0,0,1,0,0,0,1) 
//  );    

    //   Calculate the field of surface normal unit vectors.
//  surfaceVectorField surfNorm = mesh.Sf()/mesh.magSf();

    //   Calculate the LES filter width (assumes the mesh is uniform spacing throughout).
    //dimensionedScalar deltaLES 
    //(
    //     "deltaLES", 
    //	 dimensionSet(0,1,0,0,0,0,0),  
    //	 scalar(Foam::pow(mesh.V()[0],1.0/3.0))
    //);
    //deltaLES = deltaLESCoeff * deltaLES;

    //   Initialize the uStar scalar
    dimensionedScalar uStar
    (
         "uStar",
         dimensionSet(0,1,-1,0,0,0,0),
         scalar(0.0)
    );

    //   Initialize the Obukhov length
    dimensionedScalar L
    (
         "L",
         dimensionSet(0,1,0,0,0,0,0),
         scalar(1.0E30 )
    );

    //   Initialize the universal function determined from observations
    scalar phiM1 = 1.0;

    
    //   Find the distinct vertical levels of cell centers and faces normal to up
    //   for horizontal averaging
//  #include "findVerticalCellLevels.H"
//  #include "findVerticalFaceLevels.H"

    //   Find the vertical cell faces adjacent to the lower and upper boundaries
//  #include "findAdjFaces.H"

    //   Get information so that later, the driving pressure gradient can be
    //   corrected such that a specified horizontally averaged velocity at a
    //   given height is acheived.
    //#include "findWindHeight.H"

    //   Open the averaging files.
    //#include "openCellStatisticsFiles.H"
    //#include "openFaceStatisticsFiles.H"
    //#include "openABLStatisticsFiles.H"

    scalar meanTimeSum = max((runTime.value() - meanAvgStartTime), 0.0);
    scalar corrTimeSum = max((runTime.value() - corrAvgStartTime), 0.0);
    Info << "meanTimeSum = " << meanTimeSum << tab << "corrTimeSum = " << corrTimeSum << endl;



    // Begin time stepping
    Info << endl << "Entering time loop..." << endl << endl;
    label timeStep = 0;

    while (runTime.run()) 
    {
        // Read properties to see if there were changes.
        //scalar Cs(laminarTransport.lookupOrDefault<scalar>("Cs", 0.15));
        dimensionedScalar z0(laminarTransport.lookup("z0"));
        dimensionedScalar q0(laminarTransport.lookup("q0"));        
      
        // PISO and non-orthogonal correctors
        int nCorr          = mesh.solutionDict().subDict("options").lookupOrDefault<int>("nCorrectors", 3);
        int nNonOrthCorr   = mesh.solutionDict().subDict("options").lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

	// Read time step controls in system/controlDict
        #include "readTimeControls.H"
	
	// Calculate the mean and maximum CFL number in the field
        #include "CourantNo.H"
	
	// Set Delta_t if running at constant CFL
        #include "setDeltaT.H"

	// Update the simulation time (t^(n+1) = t^n + Delta_t) and print time
	runTime++;
	timeStep++;
        Info << "Time = " << runTime.timeName() << tab
	     << "Time Step = " << timeStep << endl;

	// Solve for turbulence variables
//      if(LESModel == "standardSmagorinsky")
//      {
//          #include "standardSmagorinsky.H"
//      }
        turbulence->correct();
        #include "calcStressSGS.H"


        // Update the turbine state and forces based on previous time step wind.
        if(turbineArrayOn)
        {
            turbines.update();
        }

	// Compute the Coriolis force (neglect the component in the vertical direction).
	volVectorField fCoriolis = -2.0*(Omega^U);
	forAll(fCoriolis,cellI)
	{
	     fCoriolis[cellI].z() = 0.0;
        }
	
	// Perform the predictor step.	
	fvVectorMatrix UEqn 
	(
	    fvm::ddt(U)			                   // time derivative
          + fvm::div(phi,U)		                   // convection
          + divDevR
	  + gradPd			                   // driving pressure gradient
	  //- ((rhok - 1.0)*g)                             // Boussinesq buoyancy term
	  - fCoriolis		                           // Coriolis acceleration
	);

        if (turbineArrayOn)
        {
            UEqn -= turbines.force();
        }

	UEqn.relax();

	solve
            (
	        UEqn 
	        == 
	        fvc::reconstruct
		     (
		          - fvc::snGrad(pd)*mesh.magSf()    // pressure gradient
                          +(fvc::interpolate(rhok - 1.0)*(g & mesh.Sf()))  // Boussinesq buoyancy term
			  //+(fvc::interpolate(fCoriolis) & mesh.Sf())
	             )          
            );

	
	// Perform the corrector steps.
	volScalarField rUA = 1.0/UEqn.A();
        surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
        for (int corr=0; corr<nCorr; corr++)
        {            
	    U = rUA*UEqn.H();

            surfaceScalarField phiU
            (
                (fvc::interpolate(U) & mesh.Sf())
               + fvc::ddtPhiCorr(rUA, U, phi)
            );

            // This function adjusts the boundary values of phiU such that phiU sums to zero.
            // It must be done here before the influence of gravity is added or else the sum of phiU should
            // not equal zero.
            adjustPhi(phiU, U, pd);

            //phi = phiU + (rUAf * (((fvc::interpolate(rhok) - 1.0) * (g & mesh.Sf())) + (fvc::interpolate(fCoriolis) & mesh.Sf())));
            phi = phiU + (rUAf * (((fvc::interpolate(rhok) - 1.0) * (g & mesh.Sf()))));
            //phi = (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, U, phi);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pdEqn
                (
                    //fvm::laplacian(rUA, pd) == fvc::div(phi)
                    fvm::laplacian(rUAf, pd) == fvc::div(phi)
                );

                if (pdRefOn)
                {
             //     Pout << "Pressure Ref Cell On!!!" << endl;
	            pdEqn.setReference(pdRefCell, pdRefValue);
	        }

                if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
                {
                    pdEqn.solve(mesh.solver(pd.name() + "Final"));
                }
                else
                {
                    pdEqn.solve(mesh.solver(pd.name()));
                }

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pdEqn.flux();
                }
            }

            //U -= rUA*fvc::grad(pd);
	    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
            U.correctBoundaryConditions();


	    // Calculate divergence of velocity and velocity flux and display
	    if (corr == nCorr-1)
	    {
                volScalarField divPhi = fvc::div(phi);
                volScalarField divPhiMag  = pow(pow(divPhi,2),0.5);
                scalar minLocalPhiContErr = min(divPhiMag).value();
                reduce(minLocalPhiContErr, minOp<scalar>());
                scalar maxLocalPhiContErr = max(divPhiMag).value();
                reduce(maxLocalPhiContErr, maxOp<scalar>());
                scalar avgLocalPhiContErr = divPhiMag.weightedAverage(mesh.V()).value();
                Info << "Local Flux Continuity Error:  Min " << minLocalPhiContErr << tab
                     <<                               "Max " << maxLocalPhiContErr << tab
                     <<                     "Weighted Mean " << avgLocalPhiContErr << endl;

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
        }
	
        // Correct the driving pressure gradient
        #include "correctGradPd.H"

	// Solve the temperature equation
	bool tempEqnOn = mesh.solutionDict().subDict("options").lookupOrDefault<Switch>("tempEqnOn", false);
	if (tempEqnOn)
	{
	     fvScalarMatrix TEqn
             (
                 fvm::ddt(T)
	       + fvm::div(phi,T)
               + divQ
	     );

             TEqn.solve();
	     
	     rhok = 1.0 - (T - TRef)/TRef;
	}

	// Average data if necessary.
	//if (runTime.outputTime())
	//{
	//     #include "averageFields.H"
	//}
        #include "averageFields.H"

	// Compute statistics if necessary.
	//#include "statisticsCell.H"
        //#include "statisticsFace.H"
        //#include "statisticsABL.H"

        // Write data to file if needed.
        runTime.write();
        #include "writeGradPd.H"

	// Write CPU and Clock time.
        Info << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s" << tab
             <<     "ClockTime = " << runTime.elapsedClockTime() << " s" << endl << endl;
    }

    Info << "End" << endl;

    return 0;
}

// ************************************************************************* //
