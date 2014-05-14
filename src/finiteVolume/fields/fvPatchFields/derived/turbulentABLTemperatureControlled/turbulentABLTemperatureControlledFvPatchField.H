/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::turbulentABLTemperatureControlledFvPatchField

Group
    grpGenericBoundaryConditions

Description
    Wall boundary condition for velocity at the ground in atmospheric boundary  
    layer incompressible flow in which there is wall roughness.  The velocity
    at the surface is specified such that the velocity gradient normal to the 
    surface is equal to the normal velocity gradient at the cell face opposite
    to the boundary face, but there is no wall normal flow.  The purpose of this
    is to supply a proper vertical gradient for use in the turbulence model
    equation.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        print        | print output to screen  | no          | false
    endtable

    Example of the boundary condition specification:
    \verbatim
        lower
        {
           type            turbulentABLTemperatureControlled;
           print           true;                       //-print output (true/false)
        }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    turbulentABLTemperatureControlledFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentABLTemperatureControlledFvPatchField_H
#define turbulentABLTemperatureControlledFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class turbulentABLTemperatureControlledFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class turbulentABLTemperatureControlledFvPatchField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data

        bool printOn_;
        scalar zInversion_;
        scalar widthInversion_;
        scalar Tbelow_;
        scalar Tabove_;
        scalar TgradBelow_;
        scalar TgradAbove_;
        scalar alphaTurbulent_;
        scalar fluctScale_;
        scalar zPeak_;
        label curTimeIndex_;
        label counter_;
        Field<scalar> turbField_;
        Field<scalar> turbFieldOld_;
        Random ranGen_;

public:

    //- Runtime type information
    TypeName("turbulentABLTemperatureControlled");


    // Constructors

        //- Construct from patch and internal field
        turbulentABLTemperatureControlledFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentABLTemperatureControlledFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentABLTemperatureControlledFvPatchField
        //  onto a new patch
        turbulentABLTemperatureControlledFvPatchField
        (
            const turbulentABLTemperatureControlledFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentABLTemperatureControlledFvPatchField
        (
            const turbulentABLTemperatureControlledFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentABLTemperatureControlledFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentABLTemperatureControlledFvPatchField
        (
            const turbulentABLTemperatureControlledFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentABLTemperatureControlledFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
