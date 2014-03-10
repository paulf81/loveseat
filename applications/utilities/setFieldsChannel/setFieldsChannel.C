/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    setFieldsChannel

Description
    Initializes the flow field for turbulent plane channel flow.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Read in the existing solution files.   
Info << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Get access to the input dictionary.
Info << "Reading setFieldsChannelDict dictionary" << endl;
IOdictionary setFieldsChannelDict
(
    IOobject
    (
        "setFieldsChannelDict",
        runTime.time().system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


// Read in the setFieldsABLDict entries.
scalar pi = Foam::constant::mathematical::pi;
scalar C(setFieldsChannelDict.lookupOrDefault<scalar>("C",9.0/8.0));
scalar epsilon(setFieldsChannelDict.lookupOrDefault<scalar>("epsilon",0.1));
scalar L1(setFieldsChannelDict.lookupOrDefault<scalar>("L1",4.0*pi));
scalar L2(setFieldsChannelDict.lookupOrDefault<scalar>("L2",2.0));
scalar L3(setFieldsChannelDict.lookupOrDefault<scalar>("L3",2.0*pi));
scalar Ubar(setFieldsChannelDict.lookupOrDefault<scalar>("Ubar",1.0));
word flowDirection(setFieldsChannelDict.lookupOrDefault<word>("flowDirection","x"));
word crossChannelDirection(setFieldsChannelDict.lookupOrDefault<word>("crossChannelDirection","y"));
bool updateInternalFields(setFieldsChannelDict.lookupOrDefault<bool>("updateInternalFields",true));
bool updateBoundaryFields(setFieldsChannelDict.lookupOrDefault<bool>("updateBoundaryFields",true));


// Now calculate the field quantities.

// Update the interior fields.
if (updateInternalFields)
{
    // Velocity.
    forAll(U,cellI)
    {
        scalar x1 = 0.0;
        scalar x2 = 0.0;
        scalar x3 = 0.0;
        scalar u1 = 0.0;
        scalar u2 = 0.0;
        scalar u3 = 0.0;

        if (flowDirection == "x")
        {
            x1 = mesh.C()[cellI].x();
            if (crossChannelDirection == "y")
            {
                x2 = mesh.C()[cellI].y();
                x3 = mesh.C()[cellI].z();
            }
            else if (crossChannelDirection == "z")
            {
                x2 = mesh.C()[cellI].z();
                x3 = mesh.C()[cellI].y();
            }
        }
        else if (flowDirection == "y")
        {
            x1 = mesh.C()[cellI].y();
            if (crossChannelDirection == "x")
            {
                x2 = mesh.C()[cellI].x();
                x3 = mesh.C()[cellI].z();
            }
            else if (crossChannelDirection == "z")
            {
                x2 = mesh.C()[cellI].z();
                x3 = mesh.C()[cellI].x();
            }
        }
        else if (flowDirection == "z")
        {
            x1 = mesh.C()[cellI].z();
            if (crossChannelDirection == "y")
            {
                x2 = mesh.C()[cellI].y();
                x3 = mesh.C()[cellI].x();
            }
            else if (crossChannelDirection == "x")
            {
                x2 = mesh.C()[cellI].x();
                x3 = mesh.C()[cellI].y();
            }
        }
        x1 /= L2/2.0;
        x2 /= L2/2.0;
        x3 /= L2/2.0;
        
        u1 = C * (1.0 - Foam::pow(x2,8.0)) + (epsilon * (L1/2.0) * Foam::sin(pi * x2) * Foam::cos(4.0 * pi * x1 / L1) * Foam::sin(2.0 * pi * x3 / L3));
        u2 = -epsilon * (1.0 + Foam::cos(pi * x2)) * Foam::sin(4.0 * pi * x1 / L1) * Foam::sin(2.0 * pi * x3 / L3);
        u3 = -epsilon * (L3/2.0) * Foam::sin(4.0 * pi * x1 / L1) * Foam::sin(pi * x2) * Foam::cos(2.0 * pi *x3 / L3);

        u1 *= Ubar;
        u2 *= Ubar;
        u3 *= Ubar;

        if (flowDirection == "x")
        {
            U[cellI].x() = u1;
            if (crossChannelDirection == "y")
            {
                U[cellI].y() = u2;
                U[cellI].z() = u3;
            }
            else if (crossChannelDirection == "z")
            {
                U[cellI].z() = u2;
                U[cellI].y() = u3;
            }
        }
        else if (flowDirection == "y")
        {
            U[cellI].y() = u1;
            if (crossChannelDirection == "x")
            {
                U[cellI].x() = u2;
                U[cellI].z() = u3;
            }
            else if (crossChannelDirection == "z")
            {
                U[cellI].z() = u2;
                U[cellI].x() = u3;
            }
        }
        if (flowDirection == "z")
        {
            U[cellI].z() = u1;
            if (crossChannelDirection == "y")
            {
                U[cellI].y() = u2;
                U[cellI].x() = u3;
            }
            else if (crossChannelDirection == "x")
            {
                U[cellI].x() = u2;
                U[cellI].y() = u3;
            }
        }
    }
}


// Update the boundary field.
if (updateBoundaryFields)
{
    U.correctBoundaryConditions();
}


// Write out the updated fields.
Info<< "Writing field U" << endl;
U.write();

return 0;
}


// ************************************************************************* //

