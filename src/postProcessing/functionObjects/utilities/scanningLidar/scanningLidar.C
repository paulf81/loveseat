/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "scanningLidar.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "interpolateXY.H"
//#include "fixedValueFvPatchFields.H"
//#include "zeroGradientFvPatchFields.H"
//#include "fvScalarMatrix.H"
//#include "fvmDdt.H"
//#include "fvmDiv.H"
//#include "fvcDiv.H"
//#include "fvmLaplacian.H"
//#include "fvmSup.H"
//#include "incompressible/turbulenceModel/turbulenceModel.H"
//#include "compressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(scanningLidar, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scanningLidar::scanningLidar
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    runTime_(mesh_.time()),
    active_(true),
    degRad((Foam::constant::mathematical::pi)/180.0),
    UName_("U"),
    dt(runTime_.deltaT().value()),
    time(runTime_.timeName()),
    t(runTime_.value()),
    rndGen(123456)
{
    // Read the dictionary.
    read(dict);

    // Create the sample lines.
    createBeams();

    // Get the current beam angles.
    rotationCurrent = interpolateXY(t,beamAngleTime,beamAngleRotation);
    elevationCurrent = interpolateXY(t,beamAngleTime,beamAngleElevation);

    // Rotate the beams if necessary.
    rotateLidar();

    writeVariables();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::scanningLidar::~scanningLidar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scanningLidar::read(const dictionary& dict)
{
    if (active_)
    {
        Info<< type() << ":" << nl;

        UName_ = dict.lookupOrDefault<word>("UName", "U");
        
        List<List<scalar> > beamScanPattern(dict.lookup("beamScanPattern"));
        forAll(beamScanPattern,i)
        {
           beamScanPatternTime.append(beamScanPattern[i][0]);
           beamScanPatternVector.append(vector::zero);
           beamScanPatternVector[i].x() = beamScanPattern[i][1];
           beamScanPatternVector[i].y() = beamScanPattern[i][2];
           beamScanPatternVector[i].z() = beamScanPattern[i][3];
           beamScanPatternVector[i] /= mag(beamScanPatternVector[i]);
        }

        timeBetweenScans = dict.lookupOrDefault<scalar>("timeBetweenScans", 0.0);

        beamOrigin = dict.lookupOrDefault<vector>("beamOrigin", vector::zero);
   
        beamMaxDistance =  dict.lookupOrDefault<scalar>("beamMaxDistance", 1000.0);

        beamDistribution = dict.lookup("beamDistribution");
        scalar beamDistMax = max(beamDistribution);
        forAll(beamDistribution,i)
        {
           beamDistribution[i] /= beamDistMax;
        }

        List<List<scalar> > beamAngle(dict.lookup("beamAngle"));
        
        forAll(beamAngle,i)
        {
           beamAngleTime.append(beamAngle[i][0]);
           beamAngleRotation.append(beamAngle[i][1]);
           beamAngleElevation.append(beamAngle[i][2]);
        }

        beamRotationAxis = dict.lookup("beamRotationAxis");
        beamRotationAxis /= mag(beamRotationAxis);

        beamElevationAxis = dict.lookup("beamElevationAxis");
        beamElevationAxis /= mag(beamElevationAxis);

        perturb = dict.lookupOrDefault<scalar>("perturb", 1E-4);
    }
}


void Foam::scanningLidar::createBeams()
{
    label nBeams = beamScanPatternTime.size();
    label nSamplePoints = beamDistribution.size();

    // Create the beam points.
    for(int i = 0; i < nBeams; i++)
    {
        samplePoints.append(List<vector>(nSamplePoints,vector::zero));
        for(int j = 0; j < nSamplePoints; j++)
        {
            samplePoints[i][j] = beamOrigin + beamMaxDistance*beamDistribution[j]*beamScanPatternVector[i];
        }
    }

    // Create the perturbation vectors that will be used to break ties
    // when deciding which processor a beam point lies upon.
    if (Pstream::myProcNo() == 0)
    {
        for(int i = 0; i < nBeams; i++)
        {
            perturbVectors.append(List<vector>(nSamplePoints,vector::zero));
            for(int j = 0; j < nSamplePoints; j++)
            {
                perturbVectors[i][j] = perturb*(2.0*rndGen.vector01()-vector::one);
            }
        }
    }
    Pstream::scatter(perturbVectors);

    // Build the controlCellID list of lists.
    for(int i = 0; i < nBeams; i++)
    {
        controlCellID.append(List<label>(nSamplePoints,-1));
    }
  
    // Identify the control cell IDs.
    findControlProcAndCell();
}


void Foam::scanningLidar::findControlProcAndCell()
{
    label nBeams = beamScanPatternTime.size();
    label nSamplePoints = beamDistribution.size();
    label totalSamplePoints = nBeams*nSamplePoints;
    label iter = 0;

    List<scalar> minDisLocal(totalSamplePoints,1.0E30);
    List<scalar> minDisGlobal(totalSamplePoints,1.0E30);

    for(int i = 0; i < nBeams; i++)
    {
        for(int j = 0; j < nSamplePoints; j++)
        {
            label cellID = 0;
            scalar minDis = 1.0E6;
            forAll(mesh_.C(),k)
            {
                scalar dis = mag(mesh_.C()[k] - (samplePoints[i][j] + perturbVectors[i][j]));
                if (dis <= minDis)
                {
                    cellID = k;
                    minDis = dis;
                }
            }
            minDisLocal[iter] = minDis;
            minDisGlobal[iter] = minDis;
            controlCellID[i][j] = cellID;
            iter++;
        }
    }

    Pstream::gather(minDisGlobal,minOp<List<scalar> >());
    Pstream::scatter(minDisGlobal);

    iter = 0;
    for(int i = 0; i < nBeams; i++)
    {
        for(int j = 0; j < nSamplePoints; j++)
        {
            if(minDisGlobal[iter] != minDisLocal[iter])
            {
                controlCellID[i][j] = -1;
            }
            iter++;
        }
    }
}


void Foam::scanningLidar::sampleWinds(label i)
{
}


vector Foam::scanningLidar::rotateVector(vector v, vector rotationPoint, vector axis, scalar angle)
{
    // Declare and define the rotation matrix.
    tensor RM;
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle);
    RM.xz() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract rotation point
    // off the point to be rotated.
    v = v - rotationPoint;

    // Perform the rotation.
    v = RM & v;

    // Return the rotated point to its new location relative to the rotation point.
    v = v + rotationPoint;

    return v;
}



void Foam::scanningLidar::rotateLidar()
{
    // Store the old beam angles.
    scalar rotationOld = rotationCurrent;
    scalar elevationOld = elevationCurrent;

    // Get the current beam angles.
    rotationCurrent = interpolateXY(t,beamAngleTime,beamAngleRotation);
    elevationCurrent = interpolateXY(t,beamAngleTime,beamAngleElevation);

    // Find the change in beam angle.
    scalar deltaRotation = rotationCurrent - rotationOld;
    scalar deltaElevation = elevationCurrent - elevationOld;

    // If the change in angle is finite, then perform the rotation.
    if ((mag(deltaRotation) > 0.0) || (mag(deltaElevation) > 0.0))
    {
        label nBeams = beamScanPatternTime.size();
        label nSamplePoints = beamDistribution.size();
        for(int i = 0; i < nBeams; i++)
        {
            samplePoints.append(List<vector>(nSamplePoints,vector::zero));
            for(int j = 0; j < nSamplePoints; j++)
            {
                samplePoints[i][j] = rotateVector(samplePoints[i][j],beamOrigin,beamRotationAxis,deltaRotation);
                samplePoints[i][j] = rotateVector(samplePoints[i][j],beamOrigin,beamElevationAxis,deltaElevation);
            }
        }
        beamElevationAxis = rotateVector(beamElevationAxis,vector::zero,beamRotationAxis,deltaRotation);

        // Because points moved, the owner cell, etc. must be recomputed.
        findControlProcAndCell();    
    }
}


void Foam::scanningLidar::execute()
{
    if (active_)
    {
        // Update the time information.
        time = runTime_.timeName();
        t = runTime_.value();
        dt = runTime_.deltaT().value();

        Info << type() << " output:" << endl;

        const volVectorField& U =
            mesh_.lookupObject<volVectorField>(UName_);

        Info << t << tab << dt << endl;

        Info << endl;
    }
}


void Foam::scanningLidar::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::scanningLidar::timeSet()
{
    // Do nothing
}


void Foam::scanningLidar::write()
{
    // Do nothing
}


void Foam::scanningLidar::writeVariables()
{
    Info << "name_: " << name_ << endl;
    Info << "active_: " << active_ << endl;
    Info << "UName_: " << UName_ << endl;
    Info << "dt: " << dt << endl;
    Info << "t: " << t << endl;
    Info << "beamScanPatternTime: " << beamScanPatternTime << endl;
    Info << "beamScanPatternVector: " << beamScanPatternVector << endl;
    Info << "timeBetweenScans: " << timeBetweenScans << endl;
    Info << "beamOrigin: " << beamOrigin << endl;
    Info << "beamMaxDistance: " << beamMaxDistance << endl;
    Info << "beamDistribution: " << beamDistribution << endl;
    Info << "beamAngleTime: " << beamAngleTime << endl;
    Info << "beamAngleRotation: " << beamAngleRotation << endl;
    Info << "beamAngleElevation: " << beamAngleElevation << endl;
    Info << "beamRotationAxis: " << beamRotationAxis << endl;
    Info << "beamElevationAxis: " << beamElevationAxis << endl;
    Info << "perturb: " << perturb << endl;
    Info << "perturbVectors: " << perturbVectors << endl;
    Info << "samplePoints: " << samplePoints << endl;
    Info << "sampledWindVectors: " << sampledWindVectors << endl;
    Pout << "controlCellID: " << controlCellID << endl;
    Info << "rotationCurrent: " << rotationCurrent << endl;
    Info << "elevationCurrent: " << elevationCurrent << endl;
}


// ************************************************************************* //
