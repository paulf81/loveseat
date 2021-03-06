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
    buoyantBoussinesqPimpleFoam

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        rho_{k} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "IFstream.H"
#include "OFstream.H"
#include "wallDist.H"
#include "interpolateXY.H"
#include "interpolateSplineXY.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createPostProcessingDir.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createAverageFields.H"
    #include "createSGSTurbulenceFields.H"
    #include "createGradP.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "findVerticalCellLevels.H"
    #include "findWindHeight.H"
    #include "openCellStatisticsFiles.H"
    #include "computeDivergence.H"
    #include "createDivSchemeBlendingField.H"
    //#include "openABLStatisticsFiles.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << "Starting time loop\n" << endl;

    // Update boundary conditions before starting in case anything needs
    // updating, for example after using mapFields to interpolate initial
    // field.
    U.correctBoundaryConditions();
    phi = linearInterpolate(U) & mesh.Sf();
    #include "turbulenceCorrect.H"
    T.correctBoundaryConditions();
    p_rgh.correctBoundaryConditions();


    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << tab;
        Info << "Time Step = " << runTime.timeIndex() << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        #include "updateDivSchemeBlendingField.H"

        // --- Pressure-velocity PIMPLE corrector loop
        for (pimple.start(); pimple.loop(); pimple++)
        {
            if (pimple.nOuterCorr() != 1)
            {
                p_rgh.storePrevIter();
            }

            Info << "   Predictor..." << endl;
            #include "UEqn.H"
            #include "turbulenceCorrect.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            for (int corr=0; corr<pimple.nCorr(); corr++)
            {
                Info << "   Corrector Step " << corr + 1 << "..." << endl;
                #include "pEqn.H"
                #include "turbulenceCorrect.H"
                #include "TEqn.H"
            }

            // --- Compute the velocity flux divergence
            #include "computeDivergence.H"

            // --- Update the driving pressure gradient
            #include "correctGradP.H"

            // --- Update the turbulence fields
//          if (pimple.turbCorr())
//          {
//              #include "turbulenceCorrect.H"
//              turbulence->correct();
//          }
        }   

        if (runTime.outputTime())
        {
            #include "averageFields.H"
        }

        #include "statisticsCell.H"
//      #include "statisticsFace.H"
//      #include "statisticsABL.H"

        #include "computeSGSTurbulenceFields.H"

        runTime.write();
        #include "writeGradP.H"

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End" << endl;

    return 0;
}


// ************************************************************************* //
