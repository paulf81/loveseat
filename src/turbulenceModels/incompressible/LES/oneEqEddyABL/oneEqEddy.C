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

\*---------------------------------------------------------------------------*/

#include "oneEqEddyABL.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneEqEddyABL, 0);
addToRunTimeSelectionTable(LESModel, oneEqEddyABL, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void oneEqEddyABL::updateSubGridScaleFields()
{
    nuSgs_ = ck_*sqrt(k_)*delta();
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oneEqEddyABL::oneEqEddyABL
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenEddyVisc(U, phi, transport),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),

    TName_
    (
        coeffDict_.lookupOrDefault<word>("TName","T")
    ),

    T_(U.db().lookupObject<volScalarField>(TName_)),

    g_(U.db().lookupObject<uniformDimensionedVectorField>("g"))
{
    bound(k_, kMin_);

    updateSubGridScaleFields();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oneEqEddyABL::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);

    volVectorField gradT = fvc::grad(T_);

    tmp<volScalarField> G = 2.0*nuSgs_*magSqr(symm(gradU));
    tmp<volScalarField> Q = -(1/300.0)*g_&(nuSgs_*gradT);

    tmp<fvScalarMatrix> kEqn
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(2.0*DkEff(), k_)
    ==
       G
     + Q
     - fvm::Sp(ce_*sqrt(k_)/delta(), k_)
    );

    kEqn().relax();
    kEqn().solve();

    bound(k_, kMin_);

    updateSubGridScaleFields();




  //    d/dt(k)
  //  + div(U*k)
  //  - div(2*nuSgs*grad(k))
  //    =
  //  - dev(R)*(grad(U)
  //  + (1/theta_0)*g&q
  //  - eps


}


bool oneEqEddyABL::read()
{
    if (GenEddyVisc::read())
    {
        ck_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
