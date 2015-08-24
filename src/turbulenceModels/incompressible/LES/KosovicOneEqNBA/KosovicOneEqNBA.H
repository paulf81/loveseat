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

Class
    Foam::incompressible::LESModels::KosovicOneEqNBA

Description
    Kosovic's one-equation nonlinear anisotropic backscatter SGS model.

    Kosovic, B. and Curry J., "A Large Eddy Simulation Study of a Quasi-Steady,
    Stably Stratified Atmospheric Boundary Layer," Journal of the Atmospheric
    Sciences, Vol 57, pp. 1052--1068, 2000.

SourceFiles
    KosovicOneEqNBA.C

\*---------------------------------------------------------------------------*/

#ifndef KosovicOneEqNBA_H
#define KosovicOneEqNBA_H

#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                      Class KosovicOneEqNBA Declaration
\*---------------------------------------------------------------------------*/

class KosovicOneEqNBA
:
    public LESModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cb_;
            dimensionedScalar Cs_;
            dimensionedScalar Ce_;
            dimensionedScalar Ske_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;


        // Fields

            volScalarField k_;
            volScalarField nuSgs_;
            volSymmTensorField nonlinearStress_;
            volScalarField l_;

      
        // Other model parameters or needed information.

            word TName_;
            word kappatName_;
            const volScalarField& T_;
            const uniformDimensionedVectorField& g_;
            IOdictionary transportDict_;
            dimensionedScalar TRef_;


public:

    //- Runtime type information
    TypeName("KosovicOneEqNBA");

    // Constructors

        //- Construct from components
        KosovicOneEqNBA
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~KosovicOneEqNBA()
    {}


    // Member Functions

        //- Return the SGS viscosity
        virtual tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
