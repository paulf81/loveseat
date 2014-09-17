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
    rotateFields

Description
    Rotates the fields given a rotation axis and angle.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "transformGeometricField.H"
#include "IOobjectList.H"

using namespace Foam;

template<class GeometricField>
void RotateFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const tensor& T
)
{
    // Search list of objects for volScalarFields
    IOobjectList fields(objects.lookupClass(GeometricField::typeName));

    forAllIter(IOobjectList, fields, fieldIter)
    {
        Info<< "    Rotating " << fieldIter()->name() << endl;

        GeometricField theta(*fieldIter(), mesh);
        transform(theta, dimensionedTensor(T), theta);
        theta.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    argList::validArgs.append("rotationAxis");
    argList::validArgs.append("rotationAngle");

#   include "setRootCase.H"
#   include "createTime.H"

    vector u = args.argRead<vector>(1);
    u /= mag(u);

    scalar phi = args.argRead<scalar>(2);
    phi *= Foam::constant::mathematical::pi/180.0;

    tensor T = tensor::zero;
    // from http://en.wikipedia.org/wiki/Rotation_matrix -- "Rotation matrix from axis and angle"
    T.xx() = Foam::cos(phi) + u.x()*u.x()*(1.0 - Foam::cos(phi));
    T.xy() = u.x()*u.y()*(1.0 - Foam::cos(phi)) - u.z()*Foam::sin(phi);
    T.xz() = u.x()*u.z()*(1.0 - Foam::cos(phi)) + u.y()*Foam::sin(phi);
    T.yx() = u.y()*u.x()*(1.0 - Foam::cos(phi)) + u.z()*Foam::sin(phi);
    T.yy() = Foam::cos(phi) + u.y()*u.y()*(1.0 - Foam::cos(phi));
    T.yz() = u.y()*u.z()*(1.0 - Foam::cos(phi)) - u.x()*Foam::sin(phi);
    T.zx() = u.z()*u.x()*(1.0 - Foam::cos(phi)) - u.y()*Foam::sin(phi);
    T.zy() = u.z()*u.y()*(1.0 - Foam::cos(phi)) + u.x()*Foam::sin(phi);
    T.zz() = Foam::cos(phi) + u.z()*u.z()*(1.0 - Foam::cos(phi));


//  {
//      pointIOField points
//      (
//          IOobject
//          (
//             "points",
//              runTime.findInstance(polyMesh::meshSubDir, "points"),
//              polyMesh::meshSubDir,
//              runTime,
//              IOobject::MUST_READ,
//              IOobject::NO_WRITE,
//              false
//          )
//      );

//      points = transform(T, points);

        // Set the precision of the points data to 10
//      IOstream::defaultPrecision(10);

//      Info<< "Writing points into directory " << points.path() << nl << endl;
//      points.write();
//  }


    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        RotateFields<volVectorField>(mesh, objects, T);
        RotateFields<volSphericalTensorField>(mesh, objects, T);
        RotateFields<volSymmTensorField>(mesh, objects, T);
        RotateFields<volTensorField>(mesh, objects, T);

        RotateFields<surfaceVectorField>(mesh, objects, T);
        RotateFields<surfaceSphericalTensorField>(mesh, objects, T);
        RotateFields<surfaceSymmTensorField>(mesh, objects, T);
        RotateFields<surfaceTensorField>(mesh, objects, T);
    }

    Info<< "\nEnd.\n" << endl;

    return 0;
}


// ************************************************************************* //
