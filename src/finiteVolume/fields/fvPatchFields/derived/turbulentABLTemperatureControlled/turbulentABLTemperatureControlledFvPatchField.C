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

\*---------------------------------------------------------------------------*/

#include "turbulentABLTemperatureControlledFvPatchField.H"
#include "volFields.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentABLTemperatureControlledFvPatchField::turbulentABLTemperatureControlledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    printOn_(false),
    zInversion_(700.0),
    widthInversion_(100.0),
    Tbelow_(300.0),
    Tabove_(308.0),
    TgradBelow_(0.0),
    TgradAbove_(0.003),
    alphaTurbulent_(0.1),
    fluctScale_(0.0003),
    zPeak_(100.0),
    curTimeIndex_(-1),
    counter_(0),
    turbField_(p.size()),
    turbFieldOld_(p.size()),
    ranGen_(label(0))
{}


turbulentABLTemperatureControlledFvPatchField::turbulentABLTemperatureControlledFvPatchField
(
    const turbulentABLTemperatureControlledFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    printOn_(ptf.printOn_),
    zInversion_(ptf.zInversion_),
    widthInversion_(ptf.widthInversion_),
    Tbelow_(ptf.Tbelow_),
    Tabove_(ptf.Tabove_),
    TgradBelow_(ptf.TgradBelow_),
    TgradAbove_(ptf.TgradAbove_),
    alphaTurbulent_(ptf.alphaTurbulent_),
    fluctScale_(ptf.fluctScale_),
    zPeak_(ptf.zPeak_),
    curTimeIndex_(-1),
    counter_(0),
    turbField_(ptf.turbField_,mapper),
    turbFieldOld_(ptf.turbFieldOld_,mapper),
    ranGen_(label(0))
{}


turbulentABLTemperatureControlledFvPatchField::turbulentABLTemperatureControlledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    printOn_(dict.lookupOrDefault<bool>("print", false)),
    zInversion_(readScalar(dict.lookup("zInversion"))),
    widthInversion_(readScalar(dict.lookup("widthInversion"))),
    Tbelow_(readScalar(dict.lookup("Tbelow"))),
    Tabove_(readScalar(dict.lookup("Tabove"))),
    TgradBelow_(readScalar(dict.lookup("TgradBelow"))),
    TgradAbove_(readScalar(dict.lookup("TgradAbove"))),
    alphaTurbulent_(readScalar(dict.lookup("alphaTurbulent"))),
    fluctScale_(readScalar(dict.lookup("fluctScale"))),
    zPeak_(readScalar(dict.lookup("fluctPeakZ"))),
    curTimeIndex_(-1),
    counter_(0),
    turbField_("turbField", dict, p.size()),
    turbFieldOld_(p.size()),
    ranGen_(label(0))
{}


turbulentABLTemperatureControlledFvPatchField::turbulentABLTemperatureControlledFvPatchField
(
    const turbulentABLTemperatureControlledFvPatchField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    printOn_(ptf.printOn_),
    zInversion_(ptf.zInversion_),
    widthInversion_(ptf.widthInversion_),
    Tbelow_(ptf.Tbelow_),
    Tabove_(ptf.Tabove_),
    TgradBelow_(ptf.TgradBelow_),
    TgradAbove_(ptf.TgradAbove_),
    alphaTurbulent_(ptf.alphaTurbulent_),
    fluctScale_(ptf.fluctScale_),
    zPeak_(ptf.zPeak_),
    curTimeIndex_(-1),
    counter_(0),
    turbField_(ptf.turbField_),
    turbFieldOld_(ptf.turbFieldOld_),
    ranGen_(label(0))
{}


turbulentABLTemperatureControlledFvPatchField::turbulentABLTemperatureControlledFvPatchField
(
    const turbulentABLTemperatureControlledFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    printOn_(ptf.printOn_),
    zInversion_(ptf.zInversion_),
    widthInversion_(ptf.widthInversion_),
    Tbelow_(ptf.Tbelow_),
    Tabove_(ptf.Tabove_),
    TgradBelow_(ptf.TgradBelow_),
    TgradAbove_(ptf.TgradAbove_),
    alphaTurbulent_(ptf.alphaTurbulent_),
    fluctScale_(ptf.fluctScale_),
    zPeak_(ptf.zPeak_),
    curTimeIndex_(-1),
    counter_(0),
    turbField_(ptf.turbField_),
    turbFieldOld_(ptf.turbFieldOld_),
    ranGen_(label(0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentABLTemperatureControlledFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalar pi_ = Foam::constant::mathematical::pi;



//  Compute the new inflow velocity.
    if (curTimeIndex_ != this->db().time().timeIndex())
    {


        scalarField TLocal(patch().size(),0.0);

        forAll(patch(), faceI)
        {
            scalar z = patch().Cf()[faceI].z();
            scalar T = 0.0;

            if (z < (zInversion_ - (0.5 * widthInversion_)))
            {
                T = Tbelow_ + TgradBelow_ * (z - (zInversion_ - (0.5 * widthInversion_)));
            }
            else if ( (z >= (zInversion_ - (0.5 * widthInversion_))) && 
                      (z <= (zInversion_ + (0.5 * widthInversion_))) )
            {
                scalar slope = (Tabove_ - Tbelow_) / max(1.0E-6,widthInversion_);
                T = Tbelow_ + slope * (z - (zInversion_ - (0.5 * widthInversion_)));
            }
            else if (z > (zInversion_ + (0.5 * widthInversion_)))
            {
                T = Tabove_ + TgradAbove_ * (z - (zInversion_ + (0.5 * widthInversion_)));
            }

            TLocal[faceI] = T;
        }


        Field<scalar> randomField(this->size());
        Field<scalar> mask(this->size());

        forAll(patch(), faceI)
        {
             ranGen_.randomise(randomField[faceI]);
             scalar z = patch().Cf()[faceI].z();
             mask[faceI] = 1.65 * (z/zPeak_) * Foam::exp(-0.5*Foam::pow((z/zPeak_),2));
        }

        scalar rmsCorr = sqrt(12.0*(2.0*alphaTurbulent_ - sqr(alphaTurbulent_)))/alphaTurbulent_;

        turbFieldOld_ = turbField_;

        turbField_ = (1 - alphaTurbulent_) * turbFieldOld_ + 
                     alphaTurbulent_ * mask * rmsCorr * (randomField - 0.5) * fluctScale_ * TLocal;
        TLocal += turbField_;
    
        this->operator==(TLocal);
 
        fixedValueFvPatchScalarField::updateCoeffs();
    
        curTimeIndex_ = this->db().time().timeIndex();
    }

    if (printOn_)
    {
    }   


}


void turbulentABLTemperatureControlledFvPatchField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("print")     << printOn_   << token::END_STATEMENT << nl;
    os.writeKeyword("zInversion") << zInversion_ << token::END_STATEMENT << nl;
    os.writeKeyword("widthInversion") << widthInversion_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tbelow") << Tbelow_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tabove") << Tabove_ << token::END_STATEMENT << nl;
    os.writeKeyword("TgradBelow") << TgradBelow_ << token::END_STATEMENT << nl;
    os.writeKeyword("TgradAbove") << TgradAbove_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphaTurbulent") << alphaTurbulent_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluctScale") << fluctScale_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluctPeakZ") << zPeak_ << token::END_STATEMENT << nl;
    turbField_.writeEntry("turbField", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
   fvPatchScalarField,
   turbulentABLTemperatureControlledFvPatchField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
