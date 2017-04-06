//For each turbine, yaw 1 deg/s in direction of error, if error is less than this increment do nothing

//Info << "(yaw) numTurbines = " << numTurbines << endl;
//Info << "(yaw) i = " << i << endl;
//Info << "(yaw) superInfoFromSC = " << superInfoFromSC << endl;
//Info << "(yaw) nacYaw = " << nacYaw << endl;

//Compute the yaw error
yawError = superInfoFromSC[i * superInfoLength] - nacYaw[i];


//Wrap the yaw to be +/- 180
if (yawError < -1.0 * Foam::constant::mathematical::pi )
	yawError = yawError + 2.0 * Foam::constant::mathematical::pi;
if (yawError >  1.0 * Foam::constant::mathematical::pi )
	yawError = yawError - 2.0 * Foam::constant::mathematical::pi;

// Grab the absolute value
// Note, not using abs, because it rounds to int
if (yawError >= 0.0)
	yawErrorAbs = yawError;
else
	yawErrorAbs = -1 * yawError;

// Decide if yawing, and set yaw delta at delta angle giving 1 Deg/s
if ( yawErrorAbs < Foam::constant::mathematical::pi / 180.0 * runTime_.deltaT().value() -.001) // small error
	deltaNacYaw[i] = 0.0; // Do nothing
else if  (yawError < 0)
	deltaNacYaw[i] = -1 * Foam::constant::mathematical::pi / 180.0 * runTime_.deltaT().value();  // Yaw negative
else
	deltaNacYaw[i] = Foam::constant::mathematical::pi / 180.0 * runTime_.deltaT().value();  // Yaw positive

//Info << "(deltaNacYaw) i = " << deltaNacYaw << endl;