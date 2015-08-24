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

    Cb_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb",
            coeffDict_,
            1.44
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
        symm
        (
            pow3(k_)/sqr(epsilon_)
           *(
                Ctau1_/fEta_
               *(
                    (gradU_ & gradU_)
                  + (gradU_ & gradU_)().T()
                )
              + Ctau2_/fEta_*(gradU_ & gradU_.T())
              + Ctau3_/fEta_*(gradU_.T() & gradU_)
            )
        )
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

tmp<volSymmTensorField> KosovicOneEqNBA::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)) + nonlinearStress_,
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> KosovicOneEqNBA::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_))) + nonlinearStress_
        )
    );
}


tmp<fvVectorMatrix> KosovicOneEqNBA::divDevReff(volVectorField& U) const
{
    return
    (
        fvc::div(nonlinearStress_)
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


bool KosovicOneEqNBA::read()
{
    if (RASModel::read())
    {
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        A1_.readIfPresent(coeffDict());
        A2_.readIfPresent(coeffDict());
        Ctau1_.readIfPresent(coeffDict());
        Ctau2_.readIfPresent(coeffDict());
        Ctau3_.readIfPresent(coeffDict());
        alphaKsi_.readIfPresent(coeffDict());

        kappa_.readIfPresent(coeffDict());
        E_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void KosovicOneEqNBA::correct()
{
    LESModel::correct();

    if (!turbulence_)
    {
        return;
    }

    gradU_ = fvc::grad(U_);

    // generation term
    tmp<volScalarField> S2 = symm(gradU_) && gradU_;

    volScalarField G
    (
        "RASModel::G",
        Cmu_*sqr(k_)/epsilon_*S2
      - (nonlinearStress_ && gradU_)
    );

    #include "nonLinearWallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    #include "wallDissipationI.H"

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity

    eta_ = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU_ + gradU_.T())));
    ksi_ = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU_ - gradU_.T())));
    Cmu_ = 2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_));
    fEta_ = A2_ + pow(eta_, 3.0);

    nut_ = Cmu_*sqr(k_)/epsilon_;

    #include "wallNonlinearViscosityI.H"

    nonlinearStress_ = symm
    (
        pow(k_, 3.0)/sqr(epsilon_)
       *(
            Ctau1_/fEta_
           *(
                (gradU_ & gradU_)
              + (gradU_ & gradU_)().T()
            )
          + Ctau2_/fEta_*(gradU_ & gradU_.T())
          + Ctau3_/fEta_*(gradU_.T() & gradU_)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
