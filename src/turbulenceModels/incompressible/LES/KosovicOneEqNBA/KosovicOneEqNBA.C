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
            0.36
        )
    ),

    cs_
    (
        dimensioned<scalar>
        (
           "cs",
           Foam::sqrt((8.0*(1.0 + cb_))/(27.0*Foam::sqr(Foam::constant::mathematical::pi)))
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

    leps_
    (
        IOobject
        (
            "leps",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        delta()
    ),

    ln_
    (
        IOobject
        (
            "ln",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        delta()
    ),

    ls_
    (
        IOobject
        (
            "ls",
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

void KosovicOneEqNBA::computeLengthScales()
{
    volScalarField gradTdotg = fvc::grad(T_) & g_;
    volVectorField gradUdotz = T(fvc::grad(U())) & (g_/mag(g_));
   
    Info << "gradUdotz = " << gradUdotz[3775] << endl;
    Info << "gradTdotg = " << gradTdotg[3775] << endl;
    forAll(gradTdotg,i)
    {
        // Compute buoyancy length scale.
        ln_[i] = 0.76*sqrt(k_[i])*sqrt(TRef_.value()/max(1.0E-6,mag(gradTdotg[i])));

        // Compute the shear length scale.
        ls_[i] = 2.76*sqrt(k_[i])/max(1.0E-6,sqrt(Foam::sqr(gradUdotz[i].x()) + Foam::sqr(gradUdotz[i].y())));

        // Compute the dissipation length scale.
        leps_[i] = 1.0/Foam::sqrt((1.0/(Foam::sqr(delta()[i])))+(1.0/(Foam::sqr(ln_[i])))+(1.0/(Foam::sqr(ls_[i]))));
    }
    gradTdotg.clear();
    gradUdotz.clear();
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
    // Update the molecular viscosity, and the grid-dependent length scale.
    LESModel::correct(gradU);


    // Update the stability-dependent length scale.
    KosovicOneEqNBA::computeLengthScales();


    // Use the stability-dependent and grid-dependent length scales to form the
    // turbulent Prandtl number.
    volScalarField Prt = 1.0/(1.0 + (2.0*leps_/delta()));
    Info << "Prt = " << Prt[3775] << endl;


    // Form the SGS-energy production terms, using old values of velocity and temperature.
    volSymmTensorField devB = KosovicOneEqNBA::devBeff();
    volSymmTensorField B = KosovicOneEqNBA::B();
  //tmp<volScalarField> P_shear = 2.0*nuSgs_*magSqr(symm(gradU));
    volScalarField P_shear = -(B && T(gradU));
    volScalarField P_buoyant = (1.0/TRef_)*g_&((nuSgs_/Prt)*fvc::grad(T_));
    Info << "devB = " << devB[3775] << endl;
    Info << "B = " << B[3775] << endl;
    Info << "nuSgs = " << nuSgs_[3775] << endl;
    Info << "P_shear = " << P_shear[3775] << endl;
    Info << "P_buoyant = " << P_buoyant[3775] << endl;


    // Build the SGS-energy equation matrix system.
    tmp<fvScalarMatrix> kEqn
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(2.0*DkEff(), k_)
    ==
       P_shear
     + P_buoyant
     - fvm::Sp(ce_*sqrt(k_)/leps_, k_)
    );


    // Solve the SGS-energy equation system.
    kEqn().relax();
    kEqn().solve();


    // Bound the SGS-energy to have a minimum value set by kMin_.
    bound(k_, kMin_);


    // Computes eddy viscosity.
    nuSgs_ = ce_*delta()*sqrt(k_);


    // Update the SGS thermal conductivity.
    volScalarField& kappat_ = const_cast<volScalarField&>(U().db().lookupObject<volScalarField>(kappatName_));
    kappat_ = nuSgs_/Prt;
//  kappat_.correctBoundaryConditions();   


    // Compute the nonlinear term.
    volSymmTensorField S = symm(fvc::grad(U()));
    volTensorField O = T(skew(fvc::grad(U())));
    volSymmTensorField SS = (S & S) - ((1.0/3.0) * I * (S && S));
    volTensorField SO = (S & O);
    volTensorField OS = (O & S);
    volTensorField SOmOS1 = SO - OS;
    volSymmTensorField SOmOS2 = twoSymm(S & O);
    Info << "S = " << S[3775] << endl;
    Info << "O = " << O[3775] << endl;
    Info << "SS = " << SS[3775] << endl;
    Info << "SO = " << SO[3775] << endl;
    Info << "OS = " << OS[3775] << endl;
    Info << "SOmOS1 = " << SOmOS1[3775] << endl;
    Info << "SOmOS2 = " << SOmOS2[3775] << endl;
    Info << "leps = " << leps_[3775] << endl;
    Info << "ln = " << ln_[3775] << endl;
    Info << "ls = " << ls_[3775] << endl;
    Info << "k_ = " << k_[3775] << endl;
    Info << "ce = " << ce_ << endl;
    Info << "cs = " << cs_ << endl;
    Info << "c1 = " << c1_ << endl;
    Info << "c2 = " << c2_ << endl;


    nonlinearStress_ = -ce_ * delta() * delta() * Foam::pow(cs_,(2.0/3.0)) * Foam::pow((27.0/(8.0*Foam::constant::mathematical::pi)),(1.0/3.0)) *
    (
         c1_ * ((S & S) - ((1.0/3.0) * I * (S && S)))
       + c2_ * (twoSymm(S & O))
    );
    Info << "nonlinearStress = " << nonlinearStress_[3775] << endl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
