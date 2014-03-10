        // Compute total kinetic energy.
        scalar totalKineticEnergy = 0.0;
        forAll(U,cellI)
        {
            totalKineticEnergy += 0.5 * mesh.V()[cellI] * (U[cellI] & U[cellI]);
        }
        reduce(totalKineticEnergy, sumOp<scalar>());
        Info << "Total Kinetic Energy: " << totalKineticEnergy << endl;

