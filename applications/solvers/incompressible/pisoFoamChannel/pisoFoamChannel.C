// Adapted from buoyantBoussinesqPisoFoam/pisoFoam/icoFoam by Matt Churchfield
// National Renewable Energy Laboratory
// February 27, 2012

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
    pisoFoamChannel

Description
    Transient solver for incompressible plane channel flow.
\*---------------------------------------------------------------------------*/




// Include the finite-volume CFD library.
#include "fvCFD.H"

// Include file input/output stream libraries.
#include "IFstream.H"
#include "OFstream.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




// Main program.
int main(int argc, char *argv[])
{

    // Do various initialization tasks like read mesh, create driving 
    // pressure gradient field, set time stepping, read various controls,
    // and set up for planar averaging.
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readRunControls.H"
    #include "createGradP.H"
    #include "readTimeControls.H"
    #include "readPISOControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
//  #include "findVerticalCellLevels.H"
//  #include "createCellStatistics.H"
//  #include "createAverageFields.H"


    // Begin time stepping.
    Info << endl << "Entering time loop..." << endl << endl;
    label timeStep = 0;

    while (runTime.run()) 
    {
        // Read the PISO algorithm controls in system/fvSolution.
        #include "readPISOControls.H"

	// Read time step controls in system/controlDict.
        #include "readTimeControls.H"
	
	// Calculate the mean and maximum CFL number in the field.
        #include "CourantNo.H"
	
	// Set Delta_t if running at constant CFL.
        #include "setDeltaT.H"

	// Update the simulation time (t^(n+1) = t^n + Delta_t) and print time.
	runTime++;
	timeStep++;
        Info << "Time = " << runTime.timeName() << tab
	     << "Time Step = " << timeStep << endl;


	// Perform the predictor step.
	fvVectorMatrix UEqn 
	(
	    fvm::ddt(U)
          + fvm::div(phi,U)
	  + gradP
	  - fvm::laplacian(nu,U)
	);

	solve(UEqn == -fvc::grad(p));

	
	// Perform the corrector steps.
	volScalarField rUA = 1.0/UEqn.A();
        surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

        for (int corr=0; corr<nCorr; corr++)
        {
	    U = rUA*UEqn.H();

            phi = (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, U, phi);

	    adjustPhi(phi, U, p);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi)
                );

                if (pRefOn)
                {
	            pEqn.setReference(pRefCell, pRefValue);
	        }

                if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
                {
                    pEqn.solve(mesh.solver(p.name() + "Final"));
                }
                else
                {
                    pEqn.solve(mesh.solver(p.name()));
                }

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }
	
	
        // Correct the driving pressure gradient.
        #include "correctGradP.H"

        // Calculate divergence of velocity and velocity flux and display.
        #include "computeDivergence.H"

        // Calculate the total kinetic energy.
        #include "computeTotalKineticEnergy.H"

        // Compute time-averaged statistics if necessary.
//      if (timeStatisticsOn)
//      {
//          #include "computeAverageFields.H"
//      }

	// Compute horizontal-averaged statistics if necessary.
//      if ( (horizontalStatisticsOn) && 
//           (runTime.value() >= horizontalStatisticsStartTime) && 
//           (timeStep % horizontalStatisticsFreq == 0) )
//      {
//          #include "computeCellStatistics.H"
//      }

        // Write data to file if needed.
        runTime.write();
        #include "writeGradP.H"

	// Write CPU and Clock time.
        Info << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s" << tab
             <<     "ClockTime = " << runTime.elapsedClockTime() << " s" << endl << endl;
    }

    Info << "End" << endl;

    return 0;
}

// ************************************************************************* //
