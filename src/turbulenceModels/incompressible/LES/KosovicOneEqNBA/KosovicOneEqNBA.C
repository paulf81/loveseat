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

#include "KosovicOneEqNBA.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(KosovicOneEqNBA, 0);
addToRunTimeSelectionTable(LESModel, KosovicOneEqNBA, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

KosovicOneEqNBA::KosovicOneEqNBA
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),

    cb_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cb",
            coeffDict_,
            1.44
        )
    ),

    cs_
    (
        dimensioned<scalar>
        (
           "cs",
           Foam::sqrt((8.0*(1.0 + cb_))/(27.0*Foam::constant::mathematical::pi))
        )
    ),

    ce_
    (
        dimensioned<scalar>
        (
            "ce",
            Foam::pow(8.0*Foam::constant::mathematical::pi/27.0,(1.0/3.0))*Foam::pow(cs_,(4.0/3.0))
        )
    ),

    Ske_
    (
        dimensioned<scalar>
        (
            "Ske",
            0.5
        )
    ),

    c1_
    (
        dimensioned<scalar>
        (
            "c1",
            (Foam::sqrt(960.0)*cb_)/(7.0*(1.0+cb_)*Ske_)
        )
    ),

    c2_
    (
        dimensioned<scalar>
        (
            "c2",
            -c1_
        )
    ),

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

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nonlinearStress_
    (
        "nonlinearStress",
        symm(delta()*delta()*(fvc::grad(U) & fvc::grad(U)))
    ),

    l_
    (
        IOobject
        (
            "l",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        delta()
    ),

    TName_
    (
        coeffDict_.lookupOrDefault<word>("TName","T")
    ),

    kappatName_
    (
        coeffDict_.lookupOrDefault<word>("kappatName","kappat")
    ),

    T_(U.db().lookupObject<volScalarField>(TName_)),

    g_(U.db().lookupObject<uniformDimensionedVectorField>("g")),

    transportDict_
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    TRef_(transportDict_.lookup("TRef"))

{
    bound(k_, kMin_);

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> KosovicOneEqNBA::B() const
{
    return ((2.0/3.0)*I)*k_ - nuSgs_*twoSymm(fvc::grad(U_)) + nonlinearStress_;
}


tmp<volSymmTensorField> KosovicOneEqNBA::devBeff() const
{
    return -nuEff()*dev(twoSymm(fvc::grad(U()))) + nonlinearStress_;
}


tmp<fvVectorMatrix> KosovicOneEqNBA::divDevBeff(volVectorField& U) const
{
    return
    (
        fvc::div(nonlinearStress_)
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}

void KosovicOneEqNBA::computeLengthScale()
{
    volScalarField gradTdotg = fvc::grad(T_) & g_;
    forAll(gradTdotg,i)
    {
        // neutral/unstable
        if (gradTdotg[i] >= 0.0)
        {
            l_[i] = delta()[i];
        }
        // stable
        else
        {
            l_[i] = min(delta()[i], 0.76*sqrt(k_[i])*sqrt(TRef_.value()/mag(gradTdotg[i])));
        }
    }
    gradTdotg.clear();
}


bool KosovicOneEqNBA::read()
{
    if (LESModel::read())
    {
        cb_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


void KosovicOneEqNBA::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
