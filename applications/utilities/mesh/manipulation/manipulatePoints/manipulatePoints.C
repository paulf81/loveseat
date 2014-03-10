//Adapted from transformPoints by Matt Churchfield 2/14/2012.

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
    manipulatePoints

Description
    Manipulates the mesh points in the polyMesh directory to create sinh
    mesh stretching in one direction.

Usage
    Options are:

    -distribution vector
        Applies a stretching to the grid in each direction using one of 
        the following distributions.
        -none:  keep the stretching as it is before manipulation.
        -sinh:  Apply sinh stretching.
        -tanh:  Apply tanh stretching.

        *sinh and tanh stretching assume uniform point distribution
        of input points file.

    -factor vector
        The distributions above are scaled by the given factor to increase/
        decrease grid clustering near edges of domain.

    -extents (vector vector)
        Translates and scales the domain to have the given minimum and
        maximum extents.  Input is ((xMin yMin zMin) (xMax yMax zMax)).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IStringStream.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//  Main program:

int main(int argc, char *argv[])
{

    // Get options for how to run the code.
    argList::addOption
    (
        "distribution",
        "vector",
        "stretch grid points using 'none', 'sinh', or 'tanh' distribution"
    );
    argList::addOption
    (
        "factor",
        "vector",
        "stretch factor for the stretching distribution"
    );
    argList::addOption
    (
        "extents",
        "(vectorMin vectorMax)",
        "specify the extents of the domain with '( (xMin yMin zMin) ' "
        "(xMax yMax zMax) )'"
    );


    // Find and read in points.
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    word regionName = polyMesh::defaultRegion;
    fileName meshDir;

    if (args.optionReadIfPresent("region", regionName))
    {
        meshDir = regionName/polyMesh::meshSubDir;
    }
    else
    {
        meshDir = polyMesh::meshSubDir;
    }

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshDir, "points"),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );


    // Find the global min/max x, y, and z.
    vector xyzMin = vector::zero;
    vector xyzMax = vector::zero;

    xyzMin.x() = min(points.component(vector::X));
    xyzMax.x() = max(points.component(vector::X));
    reduce(xyzMin.x(), minOp<scalar>());
    reduce(xyzMax.x(), maxOp<scalar>());

    xyzMin.y() = min(points.component(vector::Y));
    xyzMax.y() = max(points.component(vector::Y));
    reduce(xyzMin.y(), minOp<scalar>());
    reduce(xyzMax.y(), maxOp<scalar>());

    xyzMin.z() = min(points.component(vector::Z));
    xyzMax.z() = max(points.component(vector::Z));
    reduce(xyzMin.z(), minOp<scalar>());
    reduce(xyzMax.z(), maxOp<scalar>());



    // Apply stretching.
    if (args.optionFound("distribution") && args.optionFound("factor"))
    {
        List<word> distribution(args.optionLookup("distribution")());
        List<scalar> factor(args.optionLookup("factor")());
       
        // Have an error if more than three stretching types are given.
        if (distribution.size() > 3)
        {
            FatalErrorIn("manipulatePoints")
            << "Must enter no more than 3 distribution types corresponding" << endl
            << "to x, y, and z.  You entered " << distribution
            << exit(FatalError);
        }
        else if (factor.size() > 3)
        {
            FatalErrorIn("manipulatePoints")
            << "Must enter no more than 3 distribution factors corresponding" << endl
            << "to x, y, and z.  You entered " << factor
            << exit(FatalError);
        }

        else
        {
            // Apply stretching in each direction given
            forAll(distribution,i)
            {
                // Scale this coordinate direction's extents to go from -1 to 1.
                // Apply the sinh stretching
                scalar d = 2.0/(xyzMax[i] - xyzMin[i]);


                // If no stretching, do this:
                if (distribution[i] == "none")
                {
                    // Do nothing.
                }

                // If hyperbolic tangent spacing, do this:
                else if (distribution[i] == "tanh")
                {
                    if (i == 0)
                    {
                        points.replace(vector::X, (d * (points.component(vector::X) - xyzMin[i])) - 1.0);
                        points.replace(vector::X, Foam::tanh(factor[i] * points.component(vector::X)) /
                                                  Foam::tanh(factor[i])
                                      );
                    }
                    else if (i == 1)
                    {
                        points.replace(vector::Y, (d * (points.component(vector::Y) - xyzMin[i])) - 1.0);
                        points.replace(vector::Y, Foam::tanh(factor[i] * points.component(vector::Y)) /
                                                  Foam::tanh(factor[i])
                                      );
                    }
                    else if (i == 2)
                    {
                        points.replace(vector::Z, (d * (points.component(vector::Z) - xyzMin[i])) - 1.0);
                        points.replace(vector::Z, Foam::tanh(factor[i] * points.component(vector::Z)) /
                                                  Foam::tanh(factor[i])
                                      );
                    }
                }

                // If hyperbolic sine spacing, do this:
                else if (distribution[i] == "sinh")
                {
                    if (i == 0)
                    {
                        points.replace(vector::X, (d * (points.component(vector::X) - xyzMin[i])) - 1.0);
                        points.replace(vector::X, 
                                                  (0.5 * (Foam::sign(points.component(vector::X)) + 1.0) * 
                                                         (-( ( (Foam::sinh(factor[i] * (1.0 - points.component(vector::X)))) / 
                                                         (Foam::sinh(factor[i])) ) - 1.0) ) )
                                                - (0.5 * (Foam::sign(points.component(vector::X)) - 1.0) * 
                                                         ( ( ( (Foam::sinh(factor[i] * (1.0 + points.component(vector::X)))) / 
                                                         (Foam::sinh(factor[i])) ) - 1.0) ) ) 
                                      ); 
                    }
                    else if (i == 1)
                    {
                        points.replace(vector::Y, (d * (points.component(vector::Y) - xyzMin[i])) - 1.0); 
                        points.replace(vector::Y, 
                                                  (0.5 * (Foam::sign(points.component(vector::Y)) + 1.0) * 
                                                         (-( ( (Foam::sinh(factor[i] * (1.0 - points.component(vector::Y)))) / 
                                                         (Foam::sinh(factor[i])) ) - 1.0) ) )
                                                - (0.5 * (Foam::sign(points.component(vector::Y)) - 1.0) * 
                                                         ( ( ( (Foam::sinh(factor[i] * (1.0 + points.component(vector::Y)))) / 
                                                         (Foam::sinh(factor[i])) ) - 1.0) ) ) 
                                      ); 
                    }
                    else if (i == 2)
                    {
                        points.replace(vector::Z, (d * (points.component(vector::Z) - xyzMin[i])) - 1.0); 
                        points.replace(vector::Z, 
                                                  (0.5 * (Foam::sign(points.component(vector::Z)) + 1.0) * 
                                                         (-( ( (Foam::sinh(factor[i] * (1.0 - points.component(vector::Z)))) / 
                                                         (Foam::sinh(factor[i])) ) - 1.0) ) )
                                                - (0.5 * (Foam::sign(points.component(vector::Z)) - 1.0) * 
                                                         ( ( ( (Foam::sinh(factor[i] * (1.0 + points.component(vector::Z)))) / 
                                                         (Foam::sinh(factor[i])) ) - 1.0) ) ) 
                                      ); 
                    }
                }
            }
        }
    }


    // Recalculate the domain extents.
    xyzMin.x() = min(points.component(vector::X));
    xyzMax.x() = max(points.component(vector::X));
    reduce(xyzMin.x(), minOp<scalar>());
    reduce(xyzMax.x(), maxOp<scalar>());

    xyzMin.y() = min(points.component(vector::Y));
    xyzMax.y() = max(points.component(vector::Y));
    reduce(xyzMin.y(), minOp<scalar>());
    reduce(xyzMax.y(), maxOp<scalar>());

    xyzMin.z() = min(points.component(vector::Z));
    xyzMax.z() = max(points.component(vector::Z));
    reduce(xyzMin.z(), minOp<scalar>());
    reduce(xyzMax.z(), maxOp<scalar>());


    // Scale/translate to the desired domain extents.
    if (args.optionFound("extents"))
    {
        Pair<vector> extents(args.optionLookup("extents")());

        for (int i = 0; i < 3; i++)
        {
            scalar d = ((extents[1][i] - extents[0][i])/(xyzMax[i] - xyzMin[i]));
            d /= Foam::sign(d);
            if (i == 0)
            {
                points.replace(vector::X, (d * (points.component(vector::X) - xyzMin[i])) + extents[0][i] );
            }
            else if (i == 1)
            {
                points.replace(vector::Y, (d * (points.component(vector::Y) - xyzMin[i])) + extents[0][i] );
            }
            else if (i == 2)
            {
                points.replace(vector::Z, (d * (points.component(vector::Z) - xyzMin[i])) + extents[0][i] );
            }
        }
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info<< "Writing points into directory " << points.path() << nl << endl;
    points.write();

    return 0;
}


// ************************************************************************* //
