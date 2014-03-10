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
    makeFields

Description
    Given a mesh, will create initial fields.

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
Info << "Creating field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
);

Info << "Creating field T" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("T", dimensionSet(0,0,0,1,0,0,0), 0.0)
);

Info << "Creating field p_rgh" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("p_rgh", dimensionSet(0,2,-2,0,0,0,0), 0.0)
);

Info << "Creating field kappat" << endl;
volScalarField kappat
(
    IOobject
    (
        "kappat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappat", dimensionSet(0,2,-1,0,0,0,0), 0.0)
);

Info << "Creating field nuSgs" << endl;
volScalarField nuSgs
(
    IOobject
    (
        "nuSgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("nuSgs", dimensionSet(0,2,-1,0,0,0,0), 0.0)
);

Info << "Creating field qwall" << endl;
volVectorField qwall
(
    IOobject
    (
        "qwall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("qwall", dimensionSet(0,1,-1,1,0,0,0), vector::zero)
);

Info << "Creating field Rwall" << endl;
volSymmTensorField Rwall
(
    IOobject
    (
        "Rwall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("Rwall", dimensionSet(0,2,-2,0,0,0,0), symmTensor::zero)
);


// Write out the updated fields.
Info<< "Writing fields" << endl;
U.write();
T.write();
p_rgh.write();
kappat.write();
nuSgs.write();
qwall.write();
Rwall.write();

return 0;
}


// ************************************************************************* //

