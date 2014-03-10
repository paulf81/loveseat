        dimensionedVector UbarStar = U.weightedAverage(mesh.V());

	dimensionedVector pDelta = alpha * ((Ubar - UbarStar)/rUA.weightedAverage(mesh.V()) & flowDirectionVector) * flowDirectionVector;

	U += rUA*pDelta;
	phi += rUAf*(pDelta & mesh.Sf());

	gradP -= pDelta;

	Info << "Uncorrected Ubar  = "  << UbarStar.value() << tab
             << "Pressure Gradient = " << gradP.value()   << endl;
