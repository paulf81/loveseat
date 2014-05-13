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
    Foam::turbulentABLVelocityControlledFvPatchField

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
           type            turbulentABLVelocityControlled;
           print           true;                       //-print output (true/false)
        }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    turbulentABLVelocityControlledFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentABLVelocityControlledFvPatchField_H
#define turbulentABLVelocityControlledFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class turbulentABLVelocityControlledFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class turbulentABLVelocityControlledFvPatchField
:
    public fixedValueFvPatchVectorField
{

protected:

    // Protected data

        //- Print information
        bool printOn_;

        scalar zConst_;
        scalar z0_;

        vector refLoc_;
        scalar refU_;

        label avgWindow_;

        scalar alphaTurbulent_;
        scalar alphaMean_;

        vector fluctScale_;

        scalar zPeak_;
     
        scalar uStar_;
        scalar uStarOld_;
       
        scalar kappa_;
   
        scalar windDir_;
        scalar windAng_;

        //- Current time index
        label curTimeIndex_;
        label oldTimeIndex_;

        scalar minDis;
        label minDisCellID;

        label refLocControl_;
  
        label counter_;

        scalar uHistSum;

        DynamicList<scalar> uHist;

        Field<vector> turbField_;
        Field<vector> turbFieldOld_;

        Random ranGen_;



public:

    //- Runtime type information
    TypeName("turbulentABLVelocityControlled");


    // Constructors

        //- Construct from patch and internal field
        turbulentABLVelocityControlledFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentABLVelocityControlledFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentABLVelocityControlledFvPatchField
        //  onto a new patch
        turbulentABLVelocityControlledFvPatchField
        (
            const turbulentABLVelocityControlledFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentABLVelocityControlledFvPatchField
        (
            const turbulentABLVelocityControlledFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentABLVelocityControlledFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentABLVelocityControlledFvPatchField
        (
            const turbulentABLVelocityControlledFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentABLVelocityControlledFvPatchField(*this, iF)
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
