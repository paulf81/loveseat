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

\*---------------------------------------------------------------------------*/

#include "horizontalAxisWindTurbinesALM_tn.H"
#include "interpolateXY.H"

namespace Foam
{
namespace turbineModels
{

// * * * * * * * * * * * * * *  Constructor  * * * * * * * * * * * * * * * * //

horizontalAxisWindTurbinesALM_tn::horizontalAxisWindTurbinesALM_tn
(
    const volVectorField& U
)	      
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    // Set the pointer to the velocity field
    U_(U),

    // Set the degrees to radians convesion factor.
    degRad((Foam::constant::mathematical::pi)/180.0),

    // Set the revolutions/s to radians/s conversion factor.
    rpsRadSec(2.0*(Foam::constant::mathematical::pi)),

    // Set the revolutions/min to radians/s conversion factor.
    rpmRadSec(2.0*(Foam::constant::mathematical::pi)/60.0),

    // Set the time step size.
    dt(runTime_.deltaT().value()),

    // Set the current simulation time.
    time(runTime_.timeName()),
    t(runTime_.value()),

    // Set the pastFirstTimeStep flag to false
    pastFirstTimeStep(false),


    // Initialize the velocity gradient.
    gradU
    (
        IOobject
        (
            "gradU",
            time,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor("gradU",dimVelocity/dimLength,tensor::zero)
    ),

    // Initialize the body force.
    bodyForce
    (
        IOobject
        (
            "bodyForce",
            time,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce",dimForce/dimVolume/dimDensity,vector::zero)
    )




{   
    // Define dictionary that defines the turbine array.
    IOdictionary turbineArrayProperties
    (
        IOobject
        (
            "turbineArrayProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // Read in the turbine array properties dictionary.  This is the uppermost level dictionary
    // that describes where the turbines are, what kind they are, their initial state, and 
    // information about how the actuator line method is applied to each turbine.
    {
        List<word> listTemp = turbineArrayProperties.toc();
        for (int i = 0; i < listTemp.size(); i++)
        {
            if (listTemp[i] != "globalProperties")
            {
                turbineName.append(listTemp[i]);
            }
        }
    }

    numTurbines = turbineName.size();

    outputControl = turbineArrayProperties.subDict("globalProperties").lookupOrDefault<word>("outputControl","timeStep");
    outputInterval = turbineArrayProperties.subDict("globalProperties").lookupOrDefault<scalar>("outputInterval",1);
    lastOutputTime = runTime_.startTime().value();
    outputIndex = 0;

    includeNacelleSomeTrue = false;
    includeTowerSomeTrue = false;

    forAll(turbineName,i)
    {
        turbineType.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("turbineType")));
	includeNacelle.append(bool(turbineArrayProperties.subDict(turbineName[i]).lookup("includeNacelle")));
	includeTower.append(bool(turbineArrayProperties.subDict(turbineName[i]).lookup("includeTower")));
        baseLocation.append(vector(turbineArrayProperties.subDict(turbineName[i]).lookup("baseLocation")));
        numBladePoints.append(int(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("numBladePoints"))));
        numNacellePoints.append(int(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("numNacellePoints"))));
        numTowerPoints.append(int(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("numTowerPoints"))));
        bladePointDistType.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("bladePointDistType")));
        nacellePointDistType.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("nacellePointDistType")));
        towerPointDistType.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("towerPointDistType")));
        actuatorPointInterpTypeBlade.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("actuatorPointInterpTypeBlade")));
        actuatorPointInterpTypeNacelle.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("actuatorPointInterpTypeNacelle")));
        actuatorPointInterpTypeTower.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("actuatorPointInterpTypeTower")));
        actuatorUpdateType.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("actuatorUpdateType")));
        epsilonBlade.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("epsilonBlade"))));
        epsilonNacelle.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("epsilonNacelle"))));
        epsilonTower.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("epsilonTower"))));
        nacelleSampleDistance.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("nacelleSampleDistance"))));
        towerSampleDistance.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("towerSampleDistance"))));
        tipRootLossCorrType.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("tipRootLossCorrType")));
        rotationDir.append(word(turbineArrayProperties.subDict(turbineName[i]).lookup("rotationDir")));
        rotorSpeed.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("RotSpeed"))));
        rotorSpeedF.append(rotorSpeed[i]);
        speedError.append(0.0);
        intSpeedError.append(0.0);
        rotorAzimuth.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("Azimuth"))));
        generatorTorque.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("TorqueGen"))));
        bladePitch.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("Pitch"))));
        nacYaw.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("NacYaw"))));
        fluidDensity.append(scalar(readScalar(turbineArrayProperties.subDict(turbineName[i]).lookup("fluidDensity")))); 


        if(includeNacelle[i])   //  if any of the nacelles are active, set this global boolean to true.
        {
            includeNacelleSomeTrue = true;
        }
	else   // set the number of nacelle points to 1, if the nacelle is not active.
        {
            numNacellePoints[i] = 1;
        }


        if(includeTower[i])   //  if any of the nacelles are active, set this global boolean to true.
        {
            includeTowerSomeTrue = true;
        }
	else   // set the number of tower points to 1, if the tower is not active.
        {
            numTowerPoints[i] = 1;
        }
    }


    // Catalog the various types of turbines.  For example if three turbines are GE 1.5 and 
    // two turbines are Siemens 2.3 machines, then assign the GE 1.5 an ID of 0 and the Siemens
    // 2.3 an ID of 1.
    numTurbinesDistinct = 1;
    {
        turbineTypeDistinct.append(turbineType[0]);
        forAll(turbineType,i)
        {
            bool flag = false;
            for(int j = 0; j < numTurbinesDistinct; j++)
            {
                if(turbineType[i] == turbineTypeDistinct[j])
                {
                   flag = true;
                }
            }
            if(flag == false)
            {
                numTurbinesDistinct++;
                turbineTypeDistinct.append(turbineType[i]);
            }
        }
    }
    forAll(turbineType,i)
    {
        for(int j = 0; j < numTurbinesDistinct; j++)
        {
            if(turbineType[i] == turbineTypeDistinct[j])
            {
                turbineTypeID.append(j);
            }
        }
    }




    // For each distinct turbine, read in properties of that turbine from separate
    // dictionaries.

    for(int i = 0; i < numTurbinesDistinct; i++)
    {
        // Declare the turbineProperties dictionary for the ith turbine.
        IOdictionary turbineProperties
        (
            IOobject
            (
                turbineTypeDistinct[i],
                runTime_.constant(),"turbineProperties",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read in the data.
        NumBl.append(scalar(readScalar(turbineProperties.lookup("NumBl"))));
        TipRad.append(scalar(readScalar(turbineProperties.lookup("TipRad"))));
        HubRad.append(scalar(readScalar(turbineProperties.lookup("HubRad"))));
        UndSling.append(scalar(readScalar(turbineProperties.lookup("UndSling"))));
        OverHang.append(scalar(readScalar(turbineProperties.lookup("OverHang"))));
        NacelleLength.append(scalar(readScalar(turbineProperties.lookup("NacelleLength"))));
        NacelleFrontalArea.append(scalar(readScalar(turbineProperties.lookup("NacelleFrontalArea"))));
        TowerHt.append(scalar(readScalar(turbineProperties.lookup("TowerHt"))));
        Twr2Shft.append(scalar(readScalar(turbineProperties.lookup("Twr2Shft"))));
        ShftTilt.append(scalar(readScalar(turbineProperties.lookup("ShftTilt"))));
        PreCone.append(turbineProperties.lookup("PreCone"));
        GBRatio.append(scalar(readScalar(turbineProperties.lookup("GBRatio"))));
        RatedRotSpeed.append(scalar(readScalar(turbineProperties.lookup("RatedRotSpeed"))));
        GenIner.append(scalar(readScalar(turbineProperties.lookup("GenIner"))));
        HubIner.append(scalar(readScalar(turbineProperties.lookup("HubIner"))));
        BladeIner.append(scalar(readScalar(turbineProperties.lookup("BladeIner"))));
        DriveTrainIner.append(NumBl[i]*BladeIner[i] + HubIner[i] + GBRatio[i]*GBRatio[i]*GenIner[i]);
        GenTorqueControllerType.append(word(turbineProperties.lookup("GenTorqueControllerType")));
        NacYawControllerType.append(word(turbineProperties.lookup("NacYawControllerType")));
        BladePitchControllerType.append(word(turbineProperties.lookup("BladePitchControllerType")));
        RotSpeedLimiter.append(bool(readBool(turbineProperties.lookup("RotSpeedLimiter"))));
        GenTorqueRateLimiter.append(bool(readBool(turbineProperties.lookup("GenTorqueRateLimiter"))));
        NacYawRateLimiter.append(bool(readBool(turbineProperties.lookup("NacYawRateLimiter"))));
        BladePitchRateLimiter.append(bool(readBool(turbineProperties.lookup("BladePitchRateLimiter"))));
        SpeedFilterCornerFrequency.append(scalar(readScalar(turbineProperties.lookup("SpeedFilterCornerFrequency"))));

	NacelleEquivalentRadius.append(Foam::sqrt(NacelleFrontalArea[i]/Foam::constant::mathematical::pi));


        RateLimitGenTorque.append(readScalar(turbineProperties.subDict("GenTorqueControllerParams").lookup("RateLimitGenTorque")));
        if (GenTorqueControllerType[i] == "none")
        {
            // Read nothing.
        }
        else if (GenTorqueControllerType[i] == "fiveRegion")
        {
            CutInGenSpeed.append(readScalar(turbineProperties.subDict("GenTorqueControllerParams").lookup("CutInGenSpeed")));
            Region2StartGenSpeed.append(readScalar(turbineProperties.subDict("GenTorqueControllerParams").lookup("Region2StartGenSpeed")));
            Region2EndGenSpeed.append(readScalar(turbineProperties.subDict("GenTorqueControllerParams").lookup("Region2EndGenSpeed")));
            CutInGenTorque.append(readScalar(turbineProperties.subDict("GenTorqueControllerParams").lookup("CutInGenTorque")));
            RatedGenTorque.append(readScalar(turbineProperties.subDict("GenTorqueControllerParams").lookup("RatedGenTorque")));
            RateLimitGenTorque.append(readScalar(turbineProperties.subDict("GenTorqueControllerParams").lookup("RateLimitGenTorque")));
            KGen.append(readScalar(turbineProperties.subDict("GenTorqueControllerParams").lookup("KGen")));
        }
        else if (GenTorqueControllerType[i] == "speedTorqueTable")
        {
            SpeedTorqueTable.append(turbineProperties.subDict("GenTorqueControllerParams").lookup("SpeedTorqueTable"));
            DynamicList<scalar> speedInt;
            DynamicList<scalar> torqueInt;
 
            forAll(SpeedTorqueTable[i],j)
            {
                speedInt.append(SpeedTorqueTable[i][j][0]);
                torqueInt.append(SpeedTorqueTable[i][j][1]);
            }

            SpeedGenProfile.append(speedInt);
            TorqueGenProfile.append(torqueInt);

            speedInt.clear();
            torqueInt.clear();
        }
        
       

        RateLimitBladePitch.append(readScalar(turbineProperties.subDict("BladePitchControllerParams").lookup("RateLimitBladePitch")));
        if (BladePitchControllerType[i] == "none")
        {
            // Read nothing.
        }
        else if (BladePitchControllerType[i] == "PID")
        {
            PitchK.append(readScalar(turbineProperties.subDict("BladePitchControllerParams").lookup("PitchK")));
            PitchMin.append(readScalar(turbineProperties.subDict("BladePitchControllerParams").lookup("PitchMin")));
            PitchMax.append(readScalar(turbineProperties.subDict("BladePitchControllerParams").lookup("PitchMax")));
            PitchControlKP.append(readScalar(turbineProperties.subDict("BladePitchControllerParams").lookup("PitchControlKP")));
            PitchControlKI.append(readScalar(turbineProperties.subDict("BladePitchControllerParams").lookup("PitchControlKI")));
            PitchControlKD.append(readScalar(turbineProperties.subDict("BladePitchControllerParams").lookup("PitchControlKD")));
        }



        RateLimitNacYaw.append(readScalar(turbineProperties.subDict("NacYawControllerParams").lookup("RateLimitNacYaw")));
        if (NacYawControllerType[i] == "none")
        {
            // Read nothing.
        }
        else if (NacYawControllerType[i] == "timeYawTable")
        {
        }




        AirfoilType.append(turbineProperties.lookup("Airfoils"));




        BladeData.append(turbineProperties.lookup("BladeData"));
        {
           DynamicList<scalar> station;
           DynamicList<scalar> chord;
           DynamicList<scalar> twist;
           DynamicList<label> id;

           forAll(BladeData[i], j)
           {
               station.append(BladeData[i][j][0]);
               chord.append(BladeData[i][j][1]);
               twist.append(BladeData[i][j][2]);
               id.append(BladeData[i][j][3]);
           }

           BladeStation.append(station);
           BladeChord.append(chord);
           BladeTwist.append(twist);
           BladeAirfoilTypeID.append(id);

           station.clear();
           chord.clear();
           twist.clear();
           id.clear();
        }




	TowerData.append(turbineProperties.lookup("TowerData"));
        {
           DynamicList<scalar> height;
           DynamicList<scalar> chord;
           DynamicList<scalar> twist;
	   DynamicList<label> id;
           forAll(TowerData[i], j)
           {
               height.append(TowerData[i][j][0]);
               chord.append(TowerData[i][j][1]);
               twist.append(TowerData[i][j][2]);
               id.append(TowerData[i][j][3]);
           }

           TowerHeight.append(height);
           TowerChord.append(chord);
           TowerTwist.append(twist);
           TowerAirfoilTypeID.append(id);

           height.clear();
           chord.clear();
           twist.clear();
           id.clear();
        }

    }

    


    // Catalog the various distinct types of airfoils used in the various
    // distinct types of turbines.
    int numAirfoilsDistinct = 1;
    {
        airfoilTypesDistinct.append(AirfoilType[0][0]);
        forAll(AirfoilType,i)
        {
            forAll(AirfoilType[i],j)
            {
                bool flag = false;
                for(int k = 0; k < numAirfoilsDistinct; k++)
                {
                    if(AirfoilType[i][j] == airfoilTypesDistinct[k])
                    {
                        flag = true;
                    }
                }
                if(flag == false)
                {
                    numAirfoilsDistinct++;
                    airfoilTypesDistinct.append(AirfoilType[i][j]);
                }
            }
        }
    }



    // Reassign airfoil type IDs to blades of each turbine based on the global
    // distinct list of airfoils.
    forAll(BladeAirfoilTypeID,i)
    {
        forAll(BladeAirfoilTypeID[i],j)
        {
            for(int k = 0; k < numAirfoilsDistinct; k++)
            {
                if(AirfoilType[i][BladeAirfoilTypeID[i][j]] == airfoilTypesDistinct[k])
                {
                    BladeAirfoilTypeID[i][j] = k;
                    k = numAirfoilsDistinct;
                }
            }
        }
    }



    // Reassign airfoil type IDs to tower of each turbine based on the global
    // distinct list of airfoils.
    forAll(TowerAirfoilTypeID,i)
    {
        forAll(TowerAirfoilTypeID[i],j)
        {
            for(int k = 0; k < numAirfoilsDistinct; k++)
            {
                if(AirfoilType[i][TowerAirfoilTypeID[i][j]] == airfoilTypesDistinct[k])
                {
                    TowerAirfoilTypeID[i][j] = k;
                    k = numAirfoilsDistinct;
                }
            }
        }
    }
   


  

    // For each distinct airfoil, read in the bladePointLift and drag versus angle
    // of attack data.
    for(int i = 0; i < numAirfoilsDistinct; i++)
    {
        // Declare the airfoilsProperties dictionary for the ith airfoil.
        IOdictionary airfoilProperties
        (
            IOobject
            (
                airfoilTypesDistinct[i],
                runTime_.constant(),"airfoilProperties",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read in the data.
        airfoilData.append(airfoilProperties.lookup("airfoilData"));

        DynamicList<scalar> alphaInt;
        DynamicList<scalar> ClInt;
        DynamicList<scalar> CdInt;

        forAll(airfoilData[i],j)
        {
            alphaInt.append(airfoilData[i][j][0]);
            ClInt.append(airfoilData[i][j][1]);
            CdInt.append(airfoilData[i][j][2]);
        }

        airfoilAlpha.append(alphaInt);
        airfoilCl.append(ClInt);
        airfoilCd.append(CdInt);

        alphaInt.clear();
        ClInt.clear();
        CdInt.clear();
    }




    // Convert nacelle yaw from compass directions to the standard
    // convention of 0 degrees on the + x axis with positive degrees
    // in the counter-clockwise direction.
    forAll(nacYaw,i)
    {
        nacYaw[i] = compassToStandard(nacYaw[i]);
    }



    // Convert quantities in degrees into radians (dynamic lists
    // have to be done in loops).
    rotorAzimuth   = degRad * rotorAzimuth;
    rotorSpeed  = rpmRadSec * rotorSpeed;
    rotorSpeedF = rpmRadSec * rotorSpeedF;
    nacYaw    = degRad * nacYaw;
    ShftTilt  = degRad * ShftTilt;
    SpeedFilterCornerFrequency = rpsRadSec * SpeedFilterCornerFrequency;
    RatedRotSpeed = rpmRadSec * RatedRotSpeed;
    PitchK = degRad * PitchK;
    PitchMin = degRad * PitchMin;
    PitchMax = degRad * PitchMax;
    forAll(PreCone,i)
    {
        PreCone[i] = degRad * PreCone[i];
    }

    


    // Calculate tower shaft intersection and rotor apex locations. (The
    // i-index is at the turbine array level for each turbine and the j-
    // index is for each type of turbine--if all turbines are the same, j-
    // is always 0.)  The rotor apex is not yet rotated for initial yaw;
    // that is done below.
    for(int i = 0; i < numTurbines; i++)
    {
        int j = turbineTypeID[i];
        towerShaftIntersect.append(baseLocation[i]);
        towerShaftIntersect[i].z() = towerShaftIntersect[i].z() + TowerHt[j] + Twr2Shft[j];
        rotorApex.append(towerShaftIntersect[i]);
        rotorApex[i].x() = rotorApex[i].x() + ((OverHang[j] + UndSling[j]) * Foam::cos(ShftTilt[j]));
        rotorApex[i].z() = rotorApex[i].z() +  (OverHang[j] + UndSling[j]) * Foam::sin(ShftTilt[j]);
    }




    // Define the cells that can possibly be influenced by the force
    // exerted each turbine.  In otherwords, define a sphere of cell IDs
    // around each turbine that will be saved into memory so that the
    // entire domain need not be passed through when applying the force 
    // field.  (The i-index is at the turbine array level for each 
    // turbine, the j-index is for each type of turbine--if all turbines
    // are the same, j is always 0, and the k-index is at the individual
    // blade level.)
    for(int i = 0; i < numTurbines; i++)
    {
        int j = turbineTypeID[i];

        // First compute the radius of the force projection (to the radius
        // where the projection is only 0.001 its maximum value - this seems
        // recover 99.9% of the total forces when integrated).
        bladeProjectionRadius.append(epsilonBlade[i] * Foam::sqrt(Foam::log(1.0/0.001)));
        nacelleProjectionRadius.append(epsilonNacelle[i] * Foam::sqrt(Foam::log(1.0/0.001)) + NacelleEquivalentRadius[j]);
        towerProjectionRadius.append(epsilonTower[i] * Foam::sqrt(Foam::log(1.0/0.001)) + TowerChord[j][0]);

        // Calculate the sphere of influence radius (The sphere that 
	// envelops the rotor at all yaw angles).
        scalar sphereRadius = 0.0;
        forAll(PreCone[j],k)
        {
            scalar sphereRadiusI = Foam::sqrt(Foam::sqr((OverHang[j] + UndSling[j]) + TipRad[j]*Foam::sin(PreCone[j][k])) + Foam::sqr(TipRad[j]*Foam::cos(PreCone[j][k])));
            if(sphereRadiusI > sphereRadius)
            {
                sphereRadius = sphereRadiusI;
            }
        } 
        sphereRadius += max(nacelleSampleDistance[i],bladeProjectionRadius[i]);

        // Find the cells within the region of influence.
        DynamicList<label> influenceCellsI;
        forAll(U_.mesh().cells(),cellI)
        {
	    if ((includeTower[i]) && (U_.mesh().C()[cellI].z() <= towerShaftIntersect[i].z()))
            {
                if ( Foam::sqrt(Foam::sqr(U_.mesh().C()[cellI].x() - baseLocation[i].x()) + 
                                Foam::sqr(U_.mesh().C()[cellI].y() - baseLocation[i].y())) )
                {
                    influenceCellsI.append(cellI);
                }
            }
	    else if (((includeTower[i]) && (U_.mesh().C()[cellI].z() > towerShaftIntersect[i].z())) || (includeTower[i] != true))
            {
	        if (mag(U_.mesh().C()[cellI] - towerShaftIntersect[i]) <= sphereRadius)
                {
                    influenceCellsI.append(cellI);
                }
	    }
        }
        influenceCells.append(influenceCellsI);
        influenceCellsI.clear();

        // Create a list of turbines that this processor could forseeably control.
        // If influenceCells[i] is not empty, then turbine i belongs in the list.
        if (influenceCells[i].size() > 0)
        {
            turbinesControlled.append(i);
        }
    }


    
    // Create the actuator line points (not yet rotated for initial nacelle
    // yaw or initial rotor azimuth), the actuator tower points, and the
    // actuator nacelle points. i-index is at array level, j-index is
    // for the type of turbine, k-index is for each blade, and m-index is
    // for each actuator point.  Also create other important vectors, and
    // initialize the forces, blade aligned coordinate system, and
    // wind vectors to zero.
    totBladePoints = 0;
    totNacellePoints = 0;
    totTowerPoints = 0;
    for(int i = 0; i < numTurbines; i++)
    {
        int j = turbineTypeID[i];

        // Define which way the shaft points to distinguish between
        // upwind and downwind turbines.
	uvShaftDir.append(OverHang[j]/mag(OverHang[j]));

        // Define the vector along the shaft pointing in the
        // direction of the wind.
        uvShaft.append(rotorApex[i] - towerShaftIntersect[i]);
        uvShaft[i] = (uvShaft[i]/mag(uvShaft[i])) * uvShaftDir[i];

        // Define the vector aligned with the tower pointing from
        // the ground to the nacelle.
	uvTower.append(towerShaftIntersect[i] - baseLocation[i]);
	uvTower[i] = uvTower[i]/mag(uvTower[i]);

        // Now calculate the actuator section center points for each blade
        // of each turbine in the array.  All blades points will be calculated
        // at zero azimuth (blade pointing up), and then rotated to its correct
        // position before doing a global rotation to the initial azimuth of
        // the rotor.  Also calculate the radius of each point (not including coning).

        // Calculate the width of each blade actuator section.
        bladeDs.append(DynamicList<scalar>(0));
        if(bladePointDistType[i] == "uniform")
        {
            scalar actuatorWidth = (TipRad[j]-HubRad[j])/numBladePoints[i];
            for(int m = 0; m < numBladePoints[i]; m++)
            {
                bladeDs[i].append(actuatorWidth);
            }
        }
        // Add other point distribution types here, such as cosine, tanh.
	
	// Calculate the actual locations of the blade actuator sections.
        bladePoints.append(List<List<vector> >(NumBl[j], List<vector>(numBladePoints[i],vector::zero)));
        bladeSamplePoints.append(List<List<vector> >(NumBl[j], List<vector>(numBladePoints[i],vector::zero)));
        bladePointRadius.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));
        for(int k = 0; k < NumBl[j]; k++)
        {
            vector root = rotorApex[i];
            scalar beta = PreCone[j][k] - ShftTilt[j];
            root.x() = root.x() + HubRad[j]*Foam::sin(beta);
            root.z() = root.z() + HubRad[j]*Foam::cos(beta);
          //scalar dist = HubRad[j];
            scalar dist = 0.0;
            for(int m = 0; m < numBladePoints[i]; m++)
            {
               dist = dist + 0.5*bladeDs[i][m];
               bladePoints[i][k][m].x() = root.x() + dist*Foam::sin(beta);
               bladePoints[i][k][m].y() = root.y();
               bladePoints[i][k][m].z() = root.z() + dist*Foam::cos(beta);
               bladeSamplePoints[i][k][m] = bladePoints[i][k][m];
             //bladePointRadius[i][k][m] = dist;
               bladePointRadius[i][k][m] = HubRad[j] + dist;
               totBladePoints++;
               dist = dist + 0.5*bladeDs[i][m];
            }
            // Apply rotation to get blades, other than blade 1, in the right
            // place.
            if (k > 0)
            {
                for(int m = 0; m < numBladePoints[i]; m++)
                {
                    bladePoints[i][k][m] = rotatePoint(bladePoints[i][k][m], rotorApex[i], uvShaft[i], (360.0/NumBl[j])*k*degRad);
                    bladeSamplePoints[i][k][m] = rotatePoint(bladeSamplePoints[i][k][m], rotorApex[i], uvShaft[i], (360.0/NumBl[j])*k*degRad);
                }
            }
        }



        // Compute the location of the tower section center points for each
	// turbine tower in the array.

        // Calculate the width of each tower actuator section.
        towerDs.append(DynamicList<scalar>(0));
        if(towerPointDistType[i] == "uniform")
        {
            scalar actuatorWidth = TowerHt[j]/numTowerPoints[i];
            for(int m = 0; m < numTowerPoints[i]; m++)
            {
                towerDs[i].append(actuatorWidth);
            }
        }
        // Add other point distribution types here, such as cosine, tanh.

        // Compute the actual tower location points.
        towerPoints.append(List<vector>(numTowerPoints[i],vector::zero));
        towerSamplePoints.append(List<vector>(numTowerPoints[i],vector::zero));
        towerPointHeight.append(List<scalar>(numTowerPoints[i],0.0));
        {
	   towerPoints[i][0] = baseLocation[i];
	   towerSamplePoints[i][0] = towerPoints[i][0];
	   towerSamplePoints[i][0].x() -= towerSampleDistance[i];
           totTowerPoints++;
	   for(int m = 1; m < numTowerPoints[i]; m++)
           {
	      towerPoints[i][m] = baseLocation[i];
	      towerPoints[i][m].z() = towerPoints[i][m-1].z() + towerDs[i][m];
	      towerSamplePoints[i][m] = towerPoints[i][m];
	      towerSamplePoints[i][m].x() -= towerSampleDistance[i];
	      towerPointHeight[i][m] = towerPoints[i][m].z() - baseLocation[i].z();
              totTowerPoints++;
           }
        }



        // Compute the location of the nacelle section points for each turbine
        // in the array.

        // Calculate the width of each nacelle actuator section.
        nacelleDs.append(DynamicList<scalar>(0));
        if(nacellePointDistType[i] == "uniform")
        {
            scalar actuatorWidth = NacelleLength[j]/numNacellePoints[i];
            for(int m = 0; m < numNacellePoints[i]; m++)
            {
                nacelleDs[i].append(actuatorWidth);
            }
        }
        // Add other point distribution types here, such as cosine, tanh.

        // Compute the actual nacelle location points.
        nacellePoints.append(List<vector>(numNacellePoints[i],vector::zero));
        nacelleSamplePoint.append(vector::zero);
        {
	   nacellePoints[i][0] = rotorApex[i];
	   nacelleSamplePoint[i] = nacellePoints[i][0];
	   nacelleSamplePoint[i].x() -= nacelleSampleDistance[i];
           totNacellePoints++;
	   for(int m = 1; m < numNacellePoints[i]; m++)
           {
              nacellePoints[i][m] = nacellePoints[i][m-1] + nacelleDs[i][m] * uvShaft[i];
              totNacellePoints++;
           }
        }



        // Define the size of the blade, nacelle, and tower force arrays and set to zero.
        bladePointForce.append(List<List<vector> >(NumBl[j], List<vector>(numBladePoints[i],vector::zero)));
        towerPointForce.append(List<vector>(numTowerPoints[i],vector::zero));
        nacellePointForce.append(List<vector>(numNacellePoints[i],vector::zero));
  
        // Define the size of the blade, nacelle, and tower aligned vectors array and set to zero.
        bladeAlignedVectors.append(List<List<vector> >(NumBl[j],List<vector>(3,vector::zero)));

        // Define the actuator element wind vector arrays and set them to zero.
        bladeWindVectors.append(List<List<vector> >(NumBl[j],List<vector>(numBladePoints[i],vector::zero)));
        nacelleWindVector.append(vector::zero);
        towerWindVectors.append(List<vector>(numTowerPoints[i],vector::zero));

        // Define the size of the deltaNacYaw, deltaAzimuth, and deltaPitch lists and set to zero.
        deltaNacYaw.append(0.0);
        deltaAzimuth.append(0.0);

        // Define the size of the angle of attack lists and set to zero.
        bladePointAlpha.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

        // Define the size of the wind speed magnitude lists and set to zero.
        bladePointVmag.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

        // Define the size of the coefficient of bladePointLift lists and set to zero.
        bladePointCl.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

        // Define the size of the coefficient of drag lists and set to zero.
        bladePointCd.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

        // Define the size of the bladePointLift lists and set to zero.
        bladePointLift.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

        // Define the size of the drag lists and set to zero.
        drag.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

        // Define the size of the axial force lists and set to zero.
        bladePointAxialForcePoint.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

        // Define the size of the tangential force lists and set to zero.
        bladePointTangentialForcePoint.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

        // Define the size of the rotorThrust lists and set to zero.
        rotorThrust.append(0.0);

        // Define the size of the aerodynamic torque lists and set to zero.
        rotorTorque.append(0.0);

        // Define the size of the rotor power lists and set to zero.
        rotorPower.append(0.0);

        // Define the size of the cell-containing-actuator-point-sampling ID list and set to -1.
        bladeMinDisCellID.append(List<List<label> >(NumBl[j], List<label>(numBladePoints[i],-1)));
        nacelleMinDisCellID.append(-1);
        towerMinDisCellID.append(List<label>(numTowerPoints[i],-1));
    }



    // Yaw the nacelle to initial position.
    deltaNacYaw = nacYaw;
    yawNacelle();

    // Rotate the rotor to initial azimuth angle.
    deltaAzimuth =  rotorAzimuth;
    rotateBlades();  

    // Find out which processors control each actuator line point.
    findControlProcNo();

    // Compute the wind vectors at this initial time step.
    computeWindVectors();

    // Compute the blade forces due to this wind at the initial time step.
    computeBladeForce();

    // Compute the resultant body force at this initial time step.
    computeBodyForce();

    // Open the turbine data output files and print initial information.
    openOutputFiles();
    printOutputFiles();
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void horizontalAxisWindTurbinesALM_tn::rotateBlades()
{  
    // Perform rotation turbine by turbine.
    forAll(uvShaft, i)
    {	
        // Check the rotation direction first and set the local delta azimuth
	// variable accordingly.
	scalar deltaAzimuthI = 0.0;
	if (rotationDir[i] == "cw")
	{
            deltaAzimuthI =  deltaAzimuth[i];
        }
	if (rotationDir[i] == "ccw")
	{
            deltaAzimuthI = -deltaAzimuth[i];
        }

	// Rotate turbine blades, blade by blade, point by point.
	forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                bladePoints[i][j][k] = rotatePoint(bladePoints[i][j][k], rotorApex[i], uvShaft[i], deltaAzimuthI);
                bladeSamplePoints[i][j][k] = rotatePoint(bladeSamplePoints[i][j][k], rotorApex[i], uvShaft[i], deltaAzimuthI);
            }
        }   

	// Calculate the new azimuth angle and make sure it isn't
        // bigger than 2*pi.
        if (pastFirstTimeStep)
        {
	    rotorAzimuth[i] = rotorAzimuth[i] + deltaAzimuth[i];
            if (rotorAzimuth[i] >= 2.0 * Foam::constant::mathematical::pi)
            {
                rotorAzimuth[i] -= 2.0 * Foam::constant::mathematical::pi;
            }
        }
    }
}
        

void horizontalAxisWindTurbinesALM_tn::yawNacelle()
{
    // Perform rotation turbine by turbine.
    forAll(uvTower, i)
    {
	// Rotate the rotor apex first.
        rotorApex[i] = rotatePoint(rotorApex[i], towerShaftIntersect[i], uvTower[i], deltaNacYaw[i]);

	// Recompute the shaft unit vector since the shaft has rotated.
	uvShaft[i] = rotorApex[i] - towerShaftIntersect[i];
	uvShaft[i] = (uvShaft[i]/mag(uvShaft[i])) * uvShaftDir[i];
	
	// Rotate turbine blade points and velocity sampling points.
	forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                bladePoints[i][j][k] = rotatePoint(bladePoints[i][j][k], towerShaftIntersect[i], uvTower[i], deltaNacYaw[i]);
                bladeSamplePoints[i][j][k] = rotatePoint(bladeSamplePoints[i][j][k], towerShaftIntersect[i], uvTower[i], deltaNacYaw[i]);
            }
        }   

        // Rotate the nacelle points and the nacelle velocity sampling points.
        nacelleSamplePoint[i] = rotatePoint(nacelleSamplePoint[i], towerShaftIntersect[i], uvTower[i], deltaNacYaw[i]);
        forAll(nacellePoints[i], j)
        {
            nacellePoints[i][j] = rotatePoint(nacellePoints[i][j], towerShaftIntersect[i], uvTower[i], deltaNacYaw[i]);
        }

        // Rotate the tower velocity sampling points.
        forAll(towerSamplePoints[i], j)
        {
            towerSamplePoints[i][j] = rotatePoint(towerSamplePoints[i][j], towerShaftIntersect[i], uvTower[i], deltaNacYaw[i]);
        }



	// Compute the new yaw angle and make sure it isn't
        // bigger than 2*pi.
        if (pastFirstTimeStep)
        {
	    nacYaw[i] = nacYaw[i] + deltaNacYaw[i];
            if (nacYaw[i] >= 2.0 * Foam::constant::mathematical::pi)
            {
                nacYaw[i] -= 2.0 * Foam::constant::mathematical::pi;
            }
        }
    }
}


void horizontalAxisWindTurbinesALM_tn::computeRotSpeed()
{
    // Proceed turbine by turbine.
    forAll(rotorSpeed, i)
    {
        // Get the turbine type index.
        int j = turbineTypeID[i];

        // If the generator torque and blade pitch controllers are both set to "none", then
        // the rotor speed will remain fixed at its initial speed.
        if ((GenTorqueControllerType[j] == "none") && (BladePitchControllerType[j] == "none"))
        {
            // Do nothing.
        }

        // Otherwise numerically solve the dynamics of the rotor to compute the new rotor speed
        // based on the summation of aerodynamic and generator torque on the rotor.
        else
        {
            rotorSpeed[i] += (dt/DriveTrainIner[j])*(rotorTorque[i]*fluidDensity[i] - GBRatio[j]*generatorTorque[i]);
        }


        // Limit the rotor speed to be positive and such that the generator does not turn
        // faster than rated.
        if (RotSpeedLimiter[j])
        {
            # include "limiters/rotSpeedLimiter.H"
        }
 
        // Compute the change in blade azimuth angle based on the time step and current rotor speed.
        deltaAzimuth[i] = rotorSpeed[i] * dt;

    }
}


void horizontalAxisWindTurbinesALM_tn::filterRotSpeed()
{
    // Proceed turbine by turbine.
    forAll(rotorSpeedF, i)
    {
        // Get the turbine type index.
        int j = turbineTypeID[i];

        // Compute the filtering coefficient based on the corner frequency and time step.
        scalar alpha = exp(-dt * SpeedFilterCornerFrequency[j]);

        // Apply a simple recursive, single-pole, low-pass filter.
        rotorSpeedF[i] = (1.0 - alpha)*rotorSpeed[i] + alpha*rotorSpeedF[i];
    }
}


void horizontalAxisWindTurbinesALM_tn::controlGenTorque()
{
    // Proceed turbine by turbine.
    forAll(generatorTorque, i)
    {
        // Get the turbine type index.
        int j = turbineTypeID[i];

        // Get the current filtered generator speed.
        scalar genSpeedF = (rotorSpeedF[i]/rpmRadSec)*GBRatio[j];

        // Initialize the commanded generator torque variable;
        scalar generatorTorqueCommanded = generatorTorque[i];



        // Apply a controller to update the rotor speed.
	if (GenTorqueControllerType[j] == "none")
        {
            #include "controllers/genTorqueControllers/none.H"
        }

	else if (GenTorqueControllerType[j] == "fiveRegion")
        {
            #include "controllers/genTorqueControllers/fiveRegion.H"
	}

        else if (GenTorqueControllerType[j] == "speedTorqueTable")
        {
            #include "controllers/genTorqueControllers/speedTorqueTable.H"
        }

        // Limit the change in generator torque.
        if (GenTorqueRateLimiter[j])
        {
            #include "limiters/genTorqueRateLimiter.H"
        }

        // Update the pitch array.
        generatorTorque[i] = generatorTorqueCommanded;
    }
}
        

void horizontalAxisWindTurbinesALM_tn::controlNacYaw()
{
    // Proceed turbine by turbine.
    forAll(deltaNacYaw, i)
    {
        // Get the turbine type index.
        int j = turbineTypeID[i];


        
        // Apply a controller to update the nacelle yaw position.
        if (NacYawControllerType[j] == "none")
        {
            // Do nothing.
	    deltaNacYaw[i] = 0.0;
        }

        else if (NacYawControllerType[j] == "simple")
        {
            // Placeholder for when this is implemented.
        }
        
        else if (NacYawControllerType[j] == "timeYawTable")
        {
        }


        
        // Limit the change in nacelle yaw angle.
        if (NacYawRateLimiter[j])
        {
        }

    }
}
        

void horizontalAxisWindTurbinesALM_tn::controlBladePitch()
{
    // Proceed turbine by turbine.
    forAll(bladePitch, i)
    {

        // Get the turbine type index.
        int j = turbineTypeID[i];
        
        // Initialize the gain scheduling variable.
        scalar GK = 0.0;

        // Initialize the commanded pitch variable.
        scalar bladePitchCommanded = bladePitch[i]*degRad;


        // Apply a controller to update the blade pitch position.
        if (BladePitchControllerType[j] == "none")
        {
            #include "controllers/bladePitchControllers/none.H"
        }

        else if (BladePitchControllerType[j] == "PID")
        {
            #include "controllers/bladePitchControllers/PID.H"
        }

        // Apply pitch rate limiter.
        if (BladePitchRateLimiter[j])
        {
            #include "limiters/bladePitchRateLimiter.H"
        }

        // Update the pitch array.
        bladePitch[i] = bladePitchCommanded/degRad;
    }
}


void horizontalAxisWindTurbinesALM_tn::findControlProcNo()
{
    // Create a local and global list of minimum distance cells to velocity sampling 
    // points of turbines that this processor controls.  Initialize the values to huge.
    List<scalar> minDisLocalBlade(totBladePoints,1.0E30);
    List<scalar> minDisGlobalBlade(totBladePoints,1.0E30);
    List<scalar> minDisLocalNacelle(numTurbines,1.0E30);
    List<scalar> minDisGlobalNacelle(numTurbines,1.0E30);
    List<scalar> minDisLocalTower(totTowerPoints,1.0E30);
    List<scalar> minDisGlobalTower(totTowerPoints,1.0E30);


    forAll(turbinesControlled, p)
    {
        int i = turbinesControlled[p];
        int iterBlade = 0;
        int iterNacelle = 0;
        int iterTower = 0;
        if(i > 0)
        {
            for(int n = 0; n < i; n++)
            {
                iterBlade += numBladePoints[n] * NumBl[turbineTypeID[n]];
                iterNacelle += 1;
                iterTower += numTowerPoints[n];
            }
        }
        

        // Blade sampling points.
        forAll(bladeSamplePoints[i], j)
        {
            forAll(bladeSamplePoints[i][j], k)
            {
                // Find the cell that the sampling point lies within and the distance
                // from the sampling point to that cell center.
                label cellID = influenceCells[i][0];
                scalar minDis = mag(mesh_.C()[cellID] - bladeSamplePoints[i][j][k]);

                forAll(influenceCells[i], m)
                {
                    scalar dis = mag(mesh_.C()[influenceCells[i][m]] - bladeSamplePoints[i][j][k]);
                    if(dis <= minDis)
                    {
                        cellID = influenceCells[i][m];
                    }
                    minDis = mag(mesh_.C()[cellID] - bladeSamplePoints[i][j][k]);
                }
                minDisLocalBlade[iterBlade] = minDis;
                minDisGlobalBlade[iterBlade] = minDis;
                bladeMinDisCellID[i][j][k] = cellID;
                iterBlade++;
            }
        }

        // Tower sampling points.
        if(includeTowerSomeTrue)
        {
            forAll(towerSamplePoints[i],j)
            {
                label cellID = influenceCells[i][0];
                scalar minDis = mag(mesh_.C()[cellID] - towerSamplePoints[i][j]);

                forAll(influenceCells[i], m)
                {
                    scalar dis = mag(mesh_.C()[influenceCells[i][m]] - towerSamplePoints[i][j]);
                    if(dis <= minDis)
                    {
                        cellID = influenceCells[i][m];
                    }
                    minDis = mag(mesh_.C()[cellID] - towerSamplePoints[i][j]);
                }
                minDisLocalTower[iterBlade] = minDis;
                minDisGlobalTower[iterBlade] = minDis;
                towerMinDisCellID[i][j] = cellID;
                iterTower++;
            }
        }

        // Nacelle sampling point.
        if(includeNacelleSomeTrue)
        {
            label cellID = influenceCells[i][0];
            scalar minDis = mag(mesh_.C()[cellID] - nacelleSamplePoint[i]);
           
            forAll(influenceCells[i], m)
            {
                scalar dis = mag(mesh_.C()[influenceCells[i][m]] - nacelleSamplePoint[i]);
                if(dis <= minDis)
                {
                    cellID = influenceCells[i][m];
                }
                minDis = mag(mesh_.C()[cellID] - nacelleSamplePoint[i]);
            }
            minDisLocalNacelle[iterNacelle] = minDis;
            minDisGlobalNacelle[iterNacelle] = minDis;
            nacelleMinDisCellID[i] = cellID;
            iterNacelle++;
        }

    }

    // Parallel gather/scatter the global minimum distance list and reduce it by keeping 
    // only the minimum values.
    Pstream::gather(minDisGlobalBlade,minOp<List<scalar> >());
    Pstream::scatter(minDisGlobalBlade);

    if(includeNacelleSomeTrue)
    {
        Pstream::gather(minDisGlobalNacelle,minOp<List<scalar> >());
        Pstream::scatter(minDisGlobalNacelle);
    }

    if(includeTowerSomeTrue)
    {
        Pstream::gather(minDisGlobalTower,minOp<List<scalar> >());
        Pstream::scatter(minDisGlobalTower);
    }

    // Compare the global to local lists.  Where the lists agree, this processor controls
    // the actuator line point.
    forAll(turbinesControlled, p)
    {
        int i = turbinesControlled[p];
        int iterBlade = 0;
        int iterNacelle = 0;
        int iterTower = 0;
        if(i > 0)
        {
            for(int n = 0; n < i; n++)
            {
                iterBlade += numBladePoints[n] * NumBl[turbineTypeID[n]];
                iterNacelle += 1;
                iterTower += numTowerPoints[n];
            }
        }
        
        forAll(bladeSamplePoints[i], j)
        {
            forAll(bladeSamplePoints[i][j], k)
            {
                if(minDisGlobalBlade[iterBlade] != minDisLocalBlade[iterBlade])
                {
                    bladeMinDisCellID[i][j][k] = -1;
                }
                iterBlade++;
            }
        }

        if(minDisGlobalNacelle[iterNacelle] != minDisLocalNacelle[iterNacelle])
	{
            nacelleMinDisCellID[i] = -1;
	}

	forAll(towerSamplePoints[i], j)
        {
            if(minDisGlobalTower[iterTower] != minDisLocalTower[iterTower])
            {
                towerMinDisCellID[i][j] = -1;
            }
	    iterTower++;
        }


    }
}	


void horizontalAxisWindTurbinesALM_tn::computeWindVectors()
{
    // Create a list of wind velocity in x, y, z coordinates for each blade, nacelle, and tower sample point.
    List<vector> bladeWindVectorsLocal(totBladePoints,vector::zero);
    List<vector> nacelleWindVectorLocal(totNacellePoints,vector::zero);
    List<vector> towerWindVectorsLocal(totTowerPoints,vector::zero);
    



    // If linear interpolation of the velocity from the CFD mesh to the actuator
    // points is used, we need velocity gradient information.
    gradU = fvc::grad(U_);

    forAll(turbinesControlled, p)
    {
        int i = turbinesControlled[p];
        int iterBlade = 0;
	int iterNacelle = 0;
	int iterTower = 0;
        if(i > 0)
        {
            for(int n = 0; n < i; n++)
            {
                iterBlade += numBladePoints[n] * NumBl[turbineTypeID[n]];
		iterNacelle += 1;
		iterTower += numTowerPoints[n];
            }
        }
        
        forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                if(bladeMinDisCellID[i][j][k] != -1)
                {
                    vector velocity(vector::zero);
		    label cellID = bladeMinDisCellID[i][j][k];
		    vector point = bladeSamplePoints[i][j][k];

                    // If the velocity interpolation is "cellCenter", then just use 
                    // the velocity at the center of the cell within which this
                    // actuator point lies.  This is also the starting point for
		    // "linear" interpolation.  It will also be the default if an
		    // unknown interpolation keyword is given.
                    #include "velocityInterpolation/cellCenter.H"

                    // But if linear interpolation is used, add a correction based
                    // on the local velocity gradient.
                    if (actuatorPointInterpTypeBlade[i] == "linear")
                    {
                        #include "velocityInterpolation/linear.H"
                    }

		    bladeWindVectorsLocal[iterBlade] = velocity;
                }
                iterBlade++;
            }
        }

	forAll(towerPoints[i], j)
        {
            if(towerMinDisCellID[i][j] != -1)
            {
                vector velocity(vector::zero);
                label cellID = towerMinDisCellID[i][j];
	        vector point = towerSamplePoints[i][j];

                // If the velocity interpolation is "cellCenter", then just use 
                // the velocity at the center of the cell within which this
                // actuator point lies.  This is also the starting point for
	        // "linear" interpolation.  It will also be the default if an
                // unknown interpolation keyword is given.
                #include "velocityInterpolation/cellCenter.H"

                // But if linear interpolation is used, add a correction based
                // on the local velocity gradient.
                if (actuatorPointInterpTypeTower[i] == "linear")
                {
                    #include "velocityInterpolation/linear.H"
                }
                towerWindVectorsLocal[iterBlade] = velocity;
            }
	    iterTower++;
        }

	if(nacelleMinDisCellID[i] != -1)
        {
            vector velocity(vector::zero);
	    label cellID = nacelleMinDisCellID[i];
	    vector point = nacelleSamplePoint[i];

            // If the velocity interpolation is "cellCenter", then just use 
            // the velocity at the center of the cell within which this
            // actuator point lies.  This is also the starting point for
	    // "linear" interpolation.  It will also be the default if an
	    // unknown interpolation keyword is given.
            #include "velocityInterpolation/cellCenter.H"

            // But if linear interpolation is used, add a correction based
            // on the local velocity gradient.
            if (actuatorPointInterpTypeNacelle[i] == "linear")
            {
               #include "velocityInterpolation/linear.H"
            }

	    nacelleWindVectorLocal[iterBlade] = velocity;           
        }


    }

    // Perform a parallel gather of this local list to the master processor and
    // and then parallel scatter the list back out to all the processors.
    Pstream::gather(bladeWindVectorsLocal,sumOp<List<vector> >());
    Pstream::scatter(bladeWindVectorsLocal);

    if(includeNacelleSomeTrue)
    {
        Pstream::gather(nacelleWindVectorLocal,sumOp<List<vector> >());
        Pstream::scatter(nacelleWindVectorLocal);
    }

    if(includeTowerSomeTrue)
    {
        Pstream::gather(towerWindVectorsLocal,sumOp<List<vector> >());
        Pstream::scatter(towerWindVectorsLocal);
    }


    // Put the gathered/scattered wind vectors into the windVector variable.
    // Proceed turbine by turbine.
    int iterBlade = 0;
    forAll(bladeWindVectors, i)
    {
        // Proceed blade by blade.
        forAll(bladeWindVectors[i], j)
        { 
            // Proceed point by point.
            forAll(bladeWindVectors[i][j], k)
            {
                // Zero the wind vector and put in the correct velocity.
                bladeWindVectors[i][j][k] = vector::zero;
                bladeWindVectors[i][j][k] = bladeWindVectorsLocal[iterBlade];

                iterBlade++;
            }
        }
    }

    int iterNacelle = 0;
    forAll(nacelleWindVector, i)
    {
        nacelleWindVector[i] = vector::zero;
	nacelleWindVector[i] = nacelleWindVectorLocal[iterNacelle];

	iterNacelle++;
    }

    int iterTower = 0;
    forAll(towerWindVectors, i)
    {
        forAll(towerWindVectors[i], j)
        {
            towerWindVectors[i][j] = vector::zero;
	    towerWindVectors[i][j] = towerWindVectorsLocal[iterTower];

	    iterTower++;
        }
    }
}


void horizontalAxisWindTurbinesALM_tn::computeBladeForce()
{
    // The nacelle wind vector is in the Cartesian coordinate system so no need to 
    // transform it.



    // The tower nacelle wind vector is in the Cartesian coordinate system so no need
    // to transform it either.



    // Take the x,y,z wind vectors and project them into the blade coordinate system.
    // Proceed turbine by turbine.
    forAll(bladeWindVectors, i)
    {
        int n = turbineTypeID[i];

        // Proceed blade by blade.
        forAll(bladeWindVectors[i], j)
        {
            // If clockwise rotating, this vector points along the blade toward the tip.
	    // If counter-clockwise rotating, this vector points along the blade toward the root.
	    if (rotationDir[i] == "cw")
	    {
                bladeAlignedVectors[i][j][2] =   bladePoints[i][j][0] - rotorApex[i];
                bladeAlignedVectors[i][j][2] =   bladeAlignedVectors[i][j][2]/mag(bladeAlignedVectors[i][j][2]);
	    }
	    else if (rotationDir[i] == "ccw")
	    {
                bladeAlignedVectors[i][j][2] = -(bladePoints[i][j][0] - rotorApex[i]);
                bladeAlignedVectors[i][j][2] =   bladeAlignedVectors[i][j][2]/mag(bladeAlignedVectors[i][j][2]);
	    }

            // This vector points in the tangential direction opposite the turbines rotation type.  It is 
            // set up this way because it will point in the direction of oncoming flow that the blade sees 
            // due to rotation.
            bladeAlignedVectors[i][j][1] = bladeAlignedVectors[i][j][2]^uvShaft[i];
            bladeAlignedVectors[i][j][1] = bladeAlignedVectors[i][j][1]/mag(bladeAlignedVectors[i][j][1]);

            // This vector points normal to the other two and toward downwind (not exactly downwind if
            // the blade is coned).  It points in the direction of the oncoming flow due to wind that the
            // blade sees.
            bladeAlignedVectors[i][j][0] = bladeAlignedVectors[i][j][1]^bladeAlignedVectors[i][j][2];
            bladeAlignedVectors[i][j][0] = bladeAlignedVectors[i][j][0]/mag(bladeAlignedVectors[i][j][0]);
            
            // Proceed point by point.
            forAll(bladeWindVectors[i][j], k)
            {
                vector bladeWindVectorsInt = bladeWindVectors[i][j][k];

                // Zero the wind vector.
                bladeWindVectors[i][j][k] = vector::zero;

                // Now put the velocity in that cell into blade-oriented coordinates and add on the
                // velocity due to blade rotation.
                bladeWindVectors[i][j][k].x() = (bladeAlignedVectors[i][j][0] & bladeWindVectorsInt);
                bladeWindVectors[i][j][k].y() = (bladeAlignedVectors[i][j][1] & bladeWindVectorsInt) + (rotorSpeed[i] * bladePointRadius[i][j][k] * cos(PreCone[n][j]));
                bladeWindVectors[i][j][k].z() = (bladeAlignedVectors[i][j][2] & bladeWindVectorsInt);
            }
        }
    }



    // Compute the tower forces at each actuator point.
    forAll(towerWindVectors, i)
    {
        if (includeTower[i])
        {
            int m = turbineTypeID[i];

            // Set the total tower thrust of the turbine to zero.  Thrust will be summed on a tower-element-
            // wise basis.
            towerThrust[i] = 0.0;

            // Set the total tower sideways force of the turbine to zero.  The force will be summed on a tower-element-
            // wise basis.
            towerHorizontalForce[i] = 0.0;

            forAll(towerWindVectors[i], j)
            {
                // Find the flow angle.

                // Interpolate the local twist angle.
                scalar twistAng = interpolate(towerPointHeight[i][j], TowerHeight[m], TowerTwist[m]);

                // Interpolate the local chord.
                scalar chord = interpolate(towerPointHeight[i][j], TowerHeight[m], TowerChord[m]);

	        // Find the local airfoil type.
                label airfoil = interpolate(towerPointHeight[i][j], TowerHeight[m], TowerAirfoilTypeID[m]);

                // Find the local velocity magnitude composed of only the axial flow (do not include the
                // flow along the tower axis).
                towerPointVmag[i][j] = Foam::pow((Foam::pow(towerWindVectors[i][j].x(),2) + Foam::pow(towerWindVectors[i][j].y(),2)),0.5);

	        // Get the angle of the wind with respect to rotor plane tangent direction.
                scalar windAng = Foam::atan2(towerWindVectors[i][j][k].x(),towerWindVectors[i][j][k].y())/degRad;	
            }
        }
    }



    // Compute the blade forces at each actuator point.
    forAll(bladeWindVectors, i)
    {
        int m = turbineTypeID[i];

        // Set the total rotor thrust of the turbine to zero.  Thrust will be summed on a blade-element-
        // wise basis.
        rotorThrust[i] = 0.0;

        // Set the total aerodynamic torque of the turbine to zero.  Thrust will be summed on a blade-element-
        // wise basis.
        rotorTorque[i] = 0.0;

        // Proceed blade by blade.
        forAll(bladeWindVectors[i], j)
        {

            // Proceed point by point.
            forAll(bladeWindVectors[i][j], k)
            {
                // Interpolate the local twist angle.
                scalar twistAng = interpolate(bladePointRadius[i][j][k], BladeStation[m], BladeTwist[m]);

                // Interpolate the local chord.
                scalar chord = interpolate(bladePointRadius[i][j][k], BladeStation[m], BladeChord[m]);

                // Find the local airfoil type.
                label airfoil = interpolate(bladePointRadius[i][j][k], BladeStation[m], BladeAirfoilTypeID[m]);

                // Find the local velocity magnitude compose of only the axial and tangential flow (do
                // not include the radial (along blade span) flow).
                bladePointVmag[i][j][k] = Foam::pow((Foam::pow(bladeWindVectors[i][j][k].x(),2) + Foam::pow(bladeWindVectors[i][j][k].y(),2)),0.5);

                // Get the angle of the wind with respect to rotor plane tangent direction.
                scalar windAng = Foam::atan2(bladeWindVectors[i][j][k].x(),bladeWindVectors[i][j][k].y())/degRad; 

                // Angle of attack is local angle of wind with respect to rotor plane tangent minus local twist.
                bladePointAlpha[i][j][k] = windAng - twistAng - bladePitch[i];

                // Use airfoil look-up tables to get coefficient of bladePointLift and drag.
                bladePointCl[i][j][k] = interpolate(bladePointAlpha[i][j][k], airfoilAlpha[airfoil], airfoilCl[airfoil]);
                bladePointCd[i][j][k] = interpolate(bladePointAlpha[i][j][k], airfoilAlpha[airfoil], airfoilCd[airfoil]);

                // Apply tip/root-loss correction factor.
                // Tip/root-loss correction factor of Glauert.
                scalar F = 1.0;

                if(tipRootLossCorrType[i] == "none")
                {
                    F = 1.0;
                }

                else if(tipRootLossCorrType[i] == "Glauert")
                {
                    scalar g = 1.0;

                    scalar ftip  = (TipRad[m] - bladePointRadius[i][j][k])/(bladePointRadius[i][j][k] * sin(windAng*degRad));
                    scalar Ftip  = (2.0/(Foam::constant::mathematical::pi)) * acos(exp(-g * (NumBl[m] / 2.0) * ftip));

                    scalar froot = (bladePointRadius[i][j][k] - HubRad[i])/(bladePointRadius[i][j][k] * sin(windAng*degRad));
                    scalar Froot = (2.0/(Foam::constant::mathematical::pi)) * acos(exp(-g * (NumBl[m] / 2.0) * froot));

                    F = Ftip * Froot;
                }

                // Using bladePointCl, bladePointCd, wind velocity, chord, and actuator element width, calculate the
                // bladePointLift and drag per density.
                //bladePointLift[i][j][k] = 0.5 * F * bladePointCl[i][j][k] * bladePointVmag[i][j][k] * bladePointVmag[i][j][k] * chord * bladeDs[i][k];
                //drag[i][j][k] = 0.5 * F * bladePointCd[i][j][k] * bladePointVmag[i][j][k] * bladePointVmag[i][j][k] * chord * bladeDs[i][k];
                bladePointCl[i][j][k] *= F;
                bladePointCd[i][j][k] *= F;
                bladePointLift[i][j][k] = 0.5 * bladePointCl[i][j][k] * bladePointVmag[i][j][k] * bladePointVmag[i][j][k] * chord * bladeDs[i][k];
                bladePointDrag[i][j][k] = 0.5 * bladePointCd[i][j][k] * bladePointVmag[i][j][k] * bladePointVmag[i][j][k] * chord * bladeDs[i][k];

                // Make the scalar bladePointLift and drag quantities vectors in the Cartesian coordinate system.
                vector dragVector = bladeAlignedVectors[i][j][0]*bladeWindVectors[i][j][k].x() + bladeAlignedVectors[i][j][1]*bladeWindVectors[i][j][k].y();
                dragVector = dragVector/mag(dragVector);

                vector liftVector = dragVector^bladeAlignedVectors[i][j][2];
                liftVector = liftVector/mag(liftVector);

                liftVector = -bladePointLift[i][j][k] * liftVector;
                dragVector = -drag[i][j][k] * dragVector;

                // Add up bladePointLift and drag to get the resultant force/density applied to this blade element.
                bladePointForce[i][j][k] = bladePointLiftVector + dragVector;

                // Find the component of the blade element force/density in the axial (along the shaft)
                // direction.
                bladePointAxialForce[i][j][k] = -bladePointForce[i][j][k] & uvShaft[i];

                // Find the component of the blade element force/density in the tangential (torque-creating)
                // direction.
                bladePointTangentialForce[i][j][k] = bladePointForce[i][j][k] & bladeAlignedVectors[i][j][1];

                // Add this blade element's contribution to rotorThrust to the total turbine rotorThrust.
                rotorThrust[i] += bladePointAxialForce[i][j][k];

                // Add this blade element's contribution to aerodynamic torque to the total turbine aerodynamic torque.
                rotorTorque[i] += bladePointTangentialForce[i][j][k] * bladePointRadius[i][j][k] * cos(PreCone[m][j]);
            }
        }

        // Compute rotor power based on aerodynamic torque and rotation speed.
        rotorPower[i] = rotorTorque[i] * rotorSpeed[i];
    }
}


void horizontalAxisWindTurbinesALM_tn::computeBodyForce()
{  
    // Zero out the body force to begin with.
    bodyForce *= 0.0;

    // Initialize variables that are integrated forces.
    scalar rotorThrustSum = 0.0;
    scalar towerThrustSum = 0.0;
    scalar nacelleThrustSum = 0.0;

    scalar rotorTorqueSum = 0.0;

    scalar rotorThrustBodyForceSum = 0.0;
    scalar towerThrustBodyForceSum = 0.0;
    scalar nacelleThrustBodyForceSum = 0.0;

    scalar rotorTorqueBodyForceSum = 0.0;



    // Compute body force due to blades.
    forAll(bladePointForce, i)
    {
	    
        int n = turbineTypeID[i];
        
        // Proceed to compute body forces for turbine i only if there are influence cells on this processor for this turbine.
        if (influenceCells[i].size() > 0)
        {

            // Get necessary axes.
	    vector thrustVector = uvShaft[i];
	    thrustVector[2] = 0.0;
	    thrustVector = thrustVector / mag(thrustVector);
            vector verticalVector = vector::zero;
	    verticalVector[2] = 1.0;
	    vector horizontalVector = thrustVector ^ verticalVector;
	    horizontalVector = horizontalVector / mag(horizontalVector);


            // For each blade.
            forAll(bladePointForce[i], j)
            {
                // For each blade point.
                forAll(bladePointForce[i][j], k)
                {
                    // For each sphere cell.
                    forAll(influenceCells[i], m)
                    {
                        scalar dis = mag(mesh_.C()[influenceCells[i][m]] - bladePoints[i][j][k]);
                        if (dis <= bladeProjectionRadius[i])
                        {
			    scalar spreading = uniformGaussian(epsilonBlade[i], dis);
                            bodyForce[influenceCells[i][m]] += bladePointForce[i][j][k] * spreading;
                            rotorThrustBodyForceSum += (-bladePointForce[i][j][k] * spreading * mesh_.V()[influenceCells[i][m]]) & thrustVector;
                            rotorTorqueBodyForceSum += ( bladePointForce[i][j][k] * spreading * bladePointRadius[i][j][k] * cos(PreCone[n][j]) * mesh_.V()[influenceCells[i][m]]) 
				                       & bladeAlignedVectors[i][j][1];
                        }
                    }
                }  
            }
        }
        rotorThrustSum += rotorThrust[i];
        rotorTorqueSum += rotorTorque[i];
    }
    reduce(rotorThrustBodyForceSum,sumOp<scalar>());
    reduce(rotorTorqueBodyForceSum,sumOp<scalar>());

    // Print information comparing the actual rotor thrust and torque to the integrated body force.
    Info << "Rotor Thrust from Body Force = " << rotorThrustBodyForceSum << tab << "Rotor Thrust from Actuator = " << rotorThrustSum << tab
	 << "Ratio = " << rotorThrustBodyForceSum/rotorThrustSum << endl;
    Info << "Rotor Torque from Body Force = " << rotorTorqueBodyForceSum << tab << "Rotor Torque from Actuator = " << rotorTorqueSum << tab 
	 << "Ratio = " << rotorTorqueBodyForceSum/rotorTorqueSum << endl;



    // Compute body force due to tower.
    forAll(towerPointForce, i)
    {
	    
        int n = turbineTypeID[i];

        if (includeTower[i])
        {
            forAll(towerPointForce[i], j)
            {

                // Get necessary axes.
                vector thrustVector = uvShaft[i];
                thrustVector[2] = 0.0;
                thrustVector = thrustVector / mag(thrustVector);
                vector verticalVector = vector::zero;
                verticalVector[2] = 1.0;
                vector horizontalVector = thrustVector ^ verticalVector;
                horizontalVector = horizontalVector / mag(horizontalVector);

            
                forAll(influenceCells[i], m)
                {
                    scalar dis = mag(mesh_.C()[influenceCells[i][m]] - towerPoints[i][j]);
                    if (dis <= projectionRadiusTower[i])
                    {
                        scalar spreading = uniformGaussian(epsilonTower[i], dis);
			bodyForce[influenceCells[i][m]] += towerPointForce[i][j] * spreading;
			towerThrustBodyForceSum += (-towerPointForce[i][j] * spreading * mesh_.V()[influenceCells[i][m]]) & thrustVector;
                    }
                }
            }
        }
        towerThrustSum += towerThrust[i];
    }
    reduce(towerThrustBodyForceSum,sumOp<scalar>());

    // Print information comparing the actual tower thrust to the integrated body force.
    Info << "Tower Thrust from BodyForce = " << towerThrustBodyForceSum << tab << "Tower Thrust from Actuator = " << towerThrustSum << tab
	 << "Ratio = " << towerThrustBodyForceSum/towerThrustSum << endl;



    // Compute body force due to nacelle.
    forAll(nacellePointForce, i)
    {
	    
        int n = turbineTypeID[i];

        if (includeNacelle[i])
        {
            forAll(nacellePointForce[i], j)
            {

                // Get necessary axes.
                vector thrustVector = uvShaft[i];
                thrustVector[2] = 0.0;
                thrustVector = thrustVector / mag(thrustVector);
                vector verticalVector = vector::zero;
                verticalVector[2] = 1.0;
                vector horizontalVector = thrustVector ^ verticalVector;
                horizontalVector = horizontalVector / mag(horizontalVector);

            
                forAll(influenceCells[i], m)
                {
                    scalar dis = mag(mesh_.C()[influenceCells[i][m]] - nacellePoints[i][j]);
                    if (dis <= projectionRadiusTower[i])
                    {
                        scalar spreading = uniformGaussian(epsilonNacelle[i], dis);
			bodyForce[influenceCells[i][m]] += nacellePointForce[i][j] * spreading;
			nacelleThrustBodyForceSum += (-nacellePointForce[i][j] * spreading * mesh_.V()[influenceCells[i][m]]) & thrustVector;
                    }
                }
            }
        }
        nacelleThrustSum += nacelleThrust[i];
    }
    reduce(nacelleThrustBodyForceSum,sumOp<scalar>());

    // Print information comparing the actual tower thrust to the integrated body force.
    Info << "Nacelle Thrust from BodyForce = " << nacelleThrustBodyForceSum << tab << "Nacelle Thrust from Actuator = " << nacelleThrustSum << tab
	 << "Ratio = " << nacelleThrustBodyForceSum/nacelleThrustSum << endl;

}


scalar horizontalAxisWindTurbinesALM_tn::uniformGaussian(scalar epsilon, scalar d)
{
    // Compute the 3-dimensional Gaussian.
    scalar f = (1.0 / (Foam::pow(epsilon,3)*Foam::pow(Foam::constant::mathematical::pi,1.5))) * Foam::exp(-Foam::sqr(d/epsilon));
    return f;
}

scalar horizontalAxisWindTurbinesALM_tn::diskGaussian(scalar rEpsilon, scalar xEpsilon, vector u, scalar r0, vector d)
{
    // Compute a spreading function that is constant over some width radially, then dies off like a Gaussian,
    // but is also Gaussian in the axial direction.

    // Get the distances between the origin and the point in radial and axial direction.
    scalar dx = d & u;
    vector drVec = d - (dx * u);
    dx = mag(dx);
    scalar dr = mag(drVec);

    // Compute the spreading function
    scalar coeff = 1.0 / (r0 * r0 * xEpsilon * Foam::pow(Foam::constant::mathematical::pi,1.5) +
		          xEpsilon * rEpsilon * rEpsilon * Foam::pow(Foam::constant::mathematical::pi,1.5) + 
			  r0 * xEpsilon * rEpsilon * Foam::sqr(Foam::constant::mathematical::pi));
    if (dr <= r0)
    {
        scalar f = coeff * Foam::exp(-Foam::sqr(dx/xEpsilon));
    }
    else if (dr > r0)
    {
        scalar f = coeff *  Foam::exp(-Foam::sqr(dx/xEpsilon)) * Foam::exp(-Foam::sqr((dr - r0)/rEpsilon));
    }

    return f;
}


vector horizontalAxisWindTurbinesALM_tn::rotatePoint(vector point, vector rotationPoint, vector axis, scalar angle)
{
    // Declare and define the rotation matrix.
    tensor RM;
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle); 
    RM.xy() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle); 
    RM.xz() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle); 
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract rotation point
    // off the point to be rotated.
    point = point - rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation point.
    point = point + rotationPoint;

    return point;
}


scalar horizontalAxisWindTurbinesALM_tn::interpolate(scalar xNew, DynamicList<scalar>& xOld, DynamicList<scalar>& yOld)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if(diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
        }
        return yOld[indexM] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
        }
        return yOld[indexM] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}


label horizontalAxisWindTurbinesALM_tn::interpolate(scalar xNew, DynamicList<scalar>& xOld, DynamicList<label>& yOld)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if(diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
        }
        return yOld[indexM] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
        }
        return yOld[indexM] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}

scalar horizontalAxisWindTurbinesALM_tn::compassToStandard(scalar dir)
{
    dir += 180.0;
    if (dir >= 360.0)
    {
       dir -= 360.0;
    }
    dir = 90.0 - dir;
    if (dir < 0.0)
    {
        dir = dir + 360.0;
    }
    return dir;
}

scalar horizontalAxisWindTurbinesALM_tn::standardToCompass(scalar dir)
{
    dir = 90.0 - dir;
    if (dir < 0.0)
    {
        dir += 360.0;
    }
    dir += 180.0;
    if (dir >= 360.0)
    {
        dir -= 360.0;
    }
    return dir;
}
    
void horizontalAxisWindTurbinesALM_tn::update()
{
    // Update the time step size.
    dt = runTime_.deltaT().value();

    // Update the current simulation time.
    time = runTime_.timeName();
    t = runTime_.value();

    if(actuatorUpdateType[0] == "oldPosition")
    {
        // Find out which processor controls which actuator point,
        // and with that informatio sample the wind at the actuator
        // points.
        findControlProcNo();
        computeWindVectors();

        // Update the rotor state.
        filterRotSpeed();
        controlGenTorque();
        controlBladePitch();
        controlNacYaw();
        computeRotSpeed();
        rotateBlades();
        yawNacelle();
    }
    else if(actuatorUpdateType[0] == "newPosition")
    {
        // Update the rotor state.
        filterRotSpeed();
        controlGenTorque();
        controlBladePitch();
        controlNacYaw();
        computeRotSpeed();
        rotateBlades();
        yawNacelle();

        // Find out which processor controls which actuator point,
        // and with that information sample the wind at the actuator
        // points.
        findControlProcNo();
        computeWindVectors();
    }

    // Compute the blade forces.
    computeBladeForce();

    // Project the blade forces as body forces.
    computeBodyForce();

    // Print turbine output to file.
        outputIndex++;

        if (outputControl == "timeStep")
        {
            if (outputIndex >= outputInterval)
    	    {
	        outputIndex = 0;
	        printOutputFiles();
	    }
        }
        else if (outputControl == "runTime")
        {
            if ((runTime_.value() - lastOutputTime) >= outputInterval)
            {
    	        lastOutputTime += outputInterval;
	        printOutputFiles();
            }
        }
        else
        {
            printOutputFiles();
        }

    // Now that at least the first time step is finished, set pastFirstTimeStep
    // to true.
    pastFirstTimeStep = true;
}


void horizontalAxisWindTurbinesALM_tn::openOutputFiles()
{
    if (Pstream::master())
    {
        // Create the name of the root of where turbine files get ouput.
        fileName rootDir;

        if (Pstream::parRun())
        {
            rootDir = runTime_.path()/"../turbineOutput";
        }
        else
        {
            rootDir = runTime_.path()/"turbineOutput";
        }

        // Check to see if the turbineOutput directory exists; if not, create it.    
        if (!isDir(rootDir))
        {
            mkDir(rootDir);
        }

        // Check to see if the start time directory exists within the turbineOutput directory; if not, create it.  
        if (!isDir(rootDir/time))
        {
            mkDir(rootDir/time);
        }



        // Create a total aerodynamic torque file.
        //rotorTorqueFile(rootDir/time/"rotorTorque");
        rotorTorqueFile_ = new OFstream(rootDir/time/"rotorTorque");
        *rotorTorqueFile_ << "#Turbine    Time(s)    dt(s)    rotor torque (N-m)" << endl;

        // Create a generator torque file.
        generatorTorqueFile_ = new OFstream(rootDir/time/"generatorTorque");
        *generatorTorqueFile_ << "#Turbine    Time(s)    dt(s)    generator torque (N-m)" << endl;

        // Create a total rotorThrust file.
        rotorThrustFile_ = new OFstream(rootDir/time/"rotorThrust");
        *rotorThrustFile_ << "#Turbine    Time(s)    dt(s)    rotorThrust (N)" << endl;

        // Create a total power file.
        rotorPowerFile_ = new OFstream(rootDir/time/"rotorPower");
        *rotorPowerFile_ << "#Turbine    Time(s)    dt(s)    rotor power (W)" << endl;

        // Create a rotation rate file.
        rotorSpeedFile_ = new OFstream(rootDir/time/"rotorSpeed");
        *rotorSpeedFile_ << "#Turbine    Time(s)    dt(s)    rotor rotation rate(rpm)" << endl;
        
        // Create a filtered rotation rate file.
        rotorSpeedFFile_ = new OFstream(rootDir/time/"rotorSpeedFiltered");
        *rotorSpeedFFile_ << "#Turbine    Time(s)    dt(s)    filtered rotor rotation rate(rpm)" << endl;

        // Create a blade 1 azimuth angle file.
        rotorAzimuthFile_ = new OFstream(rootDir/time/"rotorAzimuth");
        *rotorAzimuthFile_ << "#Turbine    Time(s)    dt(s)    blade 1 azimuth angle (degrees)" << endl;

        // Create a blade pitch angle file.
        bladePitchFile_ = new OFstream(rootDir/time/"bladePitch");
        *bladePitchFile_ << "#Turbine    Time(s)    dt(s)    blade pitch angle (degrees)" << endl;

        // Create a nacelle yaw direction file.
        nacYawFile_ = new OFstream(rootDir/time/"nacYaw");
        *nacYawFile_ << "#Turbine    Time(s)    dt(s)    nacelle yaw angle (degrees)" << endl;

        // Create an angle of attack file.
        bladePointAlphaFile_ = new OFstream(rootDir/time/"bladePointAlpha");
        *bladePointAlphaFile_ << "#Turbine    Blade    Time(s)    dt(s)    angle-of-attack(degrees)" << endl;

        // Create a wind speed magnitude file.
        bladePointVmagFile_ = new OFstream(rootDir/time/"bladePointVmag");
        *bladePointVmagFile_ << "#Turbine    Blade    Time(s)    dt(s)    Vmag(m/s)" << endl;
    
        // Create an axial wind speed file.
        bladePointVaxialFile_ = new OFstream(rootDir/time/"bladePointVaxial");
        *bladePointVaxialFile_ << "#Turbine    Blade    Time(s)    dt(s)    Vaxial(m/s)" << endl;

        // Create a tangential wind speed file.
        bladePointVtangentialFile_ = new OFstream(rootDir/time/"bladePointVtangential");
        *bladePointVtangentialFile_ << "#Turbine    Blade    Time(s)    dt(s)    Vtangential(m/s)" << endl;

        // Create a radial wind speed file.
        bladePointVradialFile_ = new OFstream(rootDir/time/"bladePointVradial");
        *bladePointVradialFile_ << "#Turbine    Blade    Time(s)    dt(s)    Vradial(m/s)" << endl;

        // Create a coefficient of bladePointLift file.
        bladePointClFile_ = new OFstream(rootDir/time/"bladePointCl");
        *bladePointClFile_ << "#Turbine    Blade    Time(s)    dt(s)    Cl" << endl;

        // Create a coefficient of drag file.
        CdFile_ = new OFstream(rootDir/time/"bladePointCd");
        *CdFile_ << "#Turbine    Blade    Time(s)    dt(s)    Cd" << endl;

        // Create a bladePointLift file.
        bladePointLiftFile_ = new OFstream(rootDir/time/"bladePointLift");
        *bladePointLiftFile_ << "#Turbine    Blade    Time(s)    dt(s)    lift (N)" << endl;

        // Create a drag file.
        bladePointDragFile_ = new OFstream(rootDir/time/"bladePointDrag");
        *bladePointDragFile_ << "#Turbine    Blade    Time(s)    dt(s)    drag (N)" << endl;

        // Create a axial force file.
        bladePointAxialForceFile_ = new OFstream(rootDir/time/"bladePointAxialForce");
        *bladePointAxialForceFile_ << "#Turbine    Blade    Time(s)    dt(s)    axial force (N)" << endl;

        // Create a tangential force file.
        bladePointTangentialForceFile_ = new OFstream(rootDir/time/"bladePointTangentialForce");
        *bladePointTangentialForceFile_ << "#Turbine    Blade    Time(s)    dt(s)    tangential force (N)" << endl;

        // Create a x-location file.
        bladePointXFile_ = new OFstream(rootDir/time/"bladePointX");
        *bladePointXFile_ << "#Turbine    Blade    Time(s)    dt(s)    x-location(m)" << endl;

        // Create a y-location file.
        bladePointYFile_ = new OFstream(rootDir/time/"bladePointY");
        *bladePointYFile_ << "#Turbine    Blade    Time(s)    dt(s)    y-location(m)" << endl;

        // Create a z-location file.
        bladePointZFile_ = new OFstream(rootDir/time/"bladePointZ");
        *bladePointZFile_ << "#Turbine    Blade    Time(s)    dt(s)    z-location(m)" << endl;

    }
}


void horizontalAxisWindTurbinesALM_tn::printOutputFiles()
{
    if (Pstream::master())
    {
        forAll(bladePoints,i)
        {
            // Write out time and delta t.
            *rotorTorqueFile_ << i << " " << time << " " << dt << " ";
            *generatorTorqueFile_ << i << " " << time << " " << dt << " ";
            *rotorThrustFile_ << i << " " << time << " " << dt << " ";
            *rotorPowerFile_ << i << " " << time << " " << dt << " ";
            *rotorSpeedFile_ << i << " " << time << " " << dt << " ";
            *rotorSpeedFFile_ << i << " " << time << " " << dt << " ";
            *rotorAzimuthFile_ << i << " " << time << " " << dt << " ";
            *bladePitchFile_ << i << " " << time << " " << dt << " ";
            *nacYawFile_ << i << " " << time << " " << dt << " ";

            // Write out information for each turbine.
            *rotorTorqueFile_ << rotorTorque[i]*fluidDensity[i] << endl;
            *generatorTorqueFile_ << generatorTorque[i] << endl;
            *rotorThrustFile_ << rotorThrust[i]*fluidDensity[i] << endl;
            *rotorPowerFile_ << rotorPower[i]*fluidDensity[i] << endl;
            *rotorSpeedFile_ << rotorSpeed[i]/rpmRadSec << endl;
            *rotorSpeedFFile_ << rotorSpeedF[i]/rpmRadSec << endl;
            *rotorAzimuthFile_ << rotorAzimuth[i]/degRad << endl;
            *bladePitchFile_ << bladePitch[i] << endl;
            *nacYawFile_ << standardToCompass(nacYaw[i]/degRad) << endl;

            // Proceed blade by blade.
            forAll(bladePoints[i], j)
            {
                // Write out time and delta t.
                *bladePointAlphaFile_ << i << " " << j << " " << time << " " << dt << " ";
                *bladePointVmagFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointVaxialFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointVtangentialFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointVradialFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointClFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointCdFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointLiftFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointDragFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointAxialForceFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointTangentialForceFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointXFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointYFile_ << i << " " << j << " " <<  time << " " << dt << " ";
                *bladePointZFile_ << i << " " << j << " " <<  time << " " << dt << " ";

                forAll(bladePoints[i][j], k)
                {   
                    *bladePointAlphaFile_ << bladePointAlpha[i][j][k] << " ";
                    *bladePointVmagFile_ << bladePointVmag[i][j][k] << " ";
                    *bladePointVaxialFile_ << bladeWindVectors[i][j][k].x() << " ";
                    *bladePointVtangentialFile_ << bladeWindVectors[i][j][k].y() << " ";
                    *bladePointVradialFile_ << bladeWindVectors[i][j][k].z() << " ";
                    *bladePointClPointFile_ << bladePointCl[i][j][k] << " ";
                    *bladePointCdPointFile_ << bladePointCd[i][j][k] << " ";
                    *bladePointLiftFile_ << bladePointLift[i][j][k]*fluidDensity[i] << " ";
                    *bladePointDragFile_ << bladePointDrag[i][j][k]*fluidDensity[i] << " ";
                    *bladePointAxialForceFile_ << bladePointAxialForce[i][j][k]*fluidDensity[i] << " ";
                    *bladePointTangentialForceFile_ << bladePointTangentialForce[i][j][k]*fluidDensity[i] << " ";
                    *bladePointXFile_ << bladePoints[i][j][k].x() << " ";
                    *bladePointYFile_ << bladePoints[i][j][k].y() << " ";
                    *bladePointZFile_ << bladePoints[i][j][k].z() << " ";
                }
                *bladePointAlphaFile_ << endl;
                *bladePointVmagFile_ << endl;
                *bladePointVaxialFile_ << endl;
                *bladePointVtangentialFile_ << endl;
                *bladePointVradialFile_ << endl;
                *bladePointClFile_ << endl;
                *bladePointCdFile_ << endl;
                *bladePointLiftFile_ << endl;
                *bladePointDragFile_ << endl;
                *bladePointAxialForceFile_ << endl;
                *bladePointTangentialForceFile_ << endl;
                *bladePointXFile_ << endl;
                *bladePointYFile_ << endl;
                *bladePointZFile_ << endl;
            }
        }
          
        *rotorTorqueFile_ << endl;
        *generatorTorqueFile_ << endl;
        *rotorThrustFile_ << endl;
        *rotorPowerFile_ << endl;
        *rotorSpeedFile_ << endl;
        *rotorSpeedFFile_ << endl;
        *rotorAzimuthFile_ << endl;
        *bladePitchFile_ << endl;
        *nacYawFile_ << endl;

        *bladePointAlphaFile_ << endl;
        *bladePointVmagFile_ << endl;
        *bladePointVaxialFile_ << endl;
        *bladePointVtangentialFile_ << endl;
        *bladePointVradialFile_ << endl;
        *bladePointClFile_ << endl;
        *bladePointCdFile_ << endl;
        *bladePointLiftFile_ << endl;
        *bladePointDragFile_ << endl;
        *bladePointAxialForceFile_ << endl;
        *bladePointTangentialForceFile_ << endl;
        *bladePointXFile_ << endl;
        *bladePointYFile_ << endl;
        *bladePointZFile_ << endl;
    }
}
   
     
void horizontalAxisWindTurbinesALM_tn::printDebug()
{
    Info << "Print Debugging Information" << endl;
    Info << "turbineType = " << turbineType << endl;
    Info << "includeNacelle = " << includeNacelle << endl;
    Info << "includeTower = " << includeTower << endl;
    Info << "baseLocation = " << baseLocation << endl;
    Info << "numBladePoints = " << numBladePoints << endl;
    Info << "numNacellePoints = " << numNacellePoints << endl;
    Info << "numTowerPoints = " << numTowerPoints << endl;
    Info << "bladePointDistType = " << bladePointDistType << endl;
    Info << "nacellePointDistType = " << nacellePointDistType << endl;
    Info << "towerPointDistType = " << towerPointDistType << endl;
    Info << "actuatorPointInterpTypeBlade = " << actuatorPointInterpTypeBlade << endl;
    Info << "actuatorPointInterpTypeNacelle = " << actuatorPointInterpTypeNacelle << endl;
    Info << "actuatorPointInterpTypeTower = " << actuatorPointInterpTypeTower << endl;
    Info << "actuatorUpdateType = " << actuatorUpdateType << endl;
    Info << "epsilonBlade = " << epsilonBlade << endl;
    Info << "epsilonNacelle = " << epsilonNacelle << endl;
    Info << "epsilonTower = " << epsilonTower << endl;
    Info << "nacelleSampleDistance = " << nacelleSampleDistance << endl;
    Info << "towerSampleDistance = " << towerSampleDistance << endl;
    Info << "bladeProjectionRadius = " << bladeProjectionRadius << endl;
    Info << "tipRootLossCorrType = " << tipRootLossCorrType << endl;
    Info << "rotationDir = " << rotationDir << endl;
    Info << "rotorSpeed = " << rotorSpeed << endl;
    Info << "rotorSpeedF = " << rotorSpeedF << endl;
    Info << "speedError = " << speedError << endl;
    Info << "intSpeedError = " << intSpeedError << endl;
    Info << "rotorAzimuth = " << rotorAzimuth << endl;
    Info << "generatorTorque = " << generatorTorque << endl;
    Info << "bladePitch = " << bladePitch << endl;
    Info << "nacYaw = " << nacYaw << endl;
    Info << "fluidDensity = " << fluidDensity << endl << endl << endl;
    
    Info << "numTurbinesDistinct = " << numTurbinesDistinct << endl;
    Info << "turbineTypeDistinct = " << turbineTypeDistinct << endl;
    Info << "turbineTypeID = " << turbineTypeID << endl << endl << endl;;

    Info << "NumBl = " << NumBl << endl;
    Info << "TipRad = " << TipRad << endl;
    Info << "HubRad = " << HubRad << endl;
    Info << "UndSling = " << UndSling << endl;
    Info << "OverHang = " << OverHang << endl;
    Info << "NacelleLength = " << NacelleLength << endl;
    Info << "TowerHt = " << TowerHt << endl;
    Info << "Twr2Shft = " << Twr2Shft << endl;
    Info << "ShftTilt = " << ShftTilt << endl;
    Info << "PreCone = " << PreCone << endl;
    Info << "GBRatio = " << GBRatio << endl;
    Info << "RatedRotSpeed = " << RatedRotSpeed << endl;
    Info << "HubIner = " << HubIner << endl;
    Info << "GenIner = " << GenIner << endl;
    Info << "BladeIner = " << BladeIner << endl;
    Info << "GenTorqueControllerType = " << GenTorqueControllerType << endl;
    Info << "NacYawControllerType = " << NacYawControllerType << endl;
    Info << "BladePitchControllerType = " << BladePitchControllerType << endl;
    Info << "RotSpeedLimiter = " << RotSpeedLimiter << endl;
    Info << "GenTorqueRateLimiter = " << GenTorqueRateLimiter << endl;
    Info << "NacYawRateLimiter = " << NacYawRateLimiter << endl;
    Info << "BladePitchRateLimiter = " << BladePitchRateLimiter << endl;
    Info << "SpeedFilterCornerFrequency = " << SpeedFilterCornerFrequency << endl;
    Info << "AirfoilType = " << AirfoilType << endl;
    Info << "BladeData = " << BladeData << endl;
    Info << "BladeStation = " << BladeStation << endl;
    Info << "BladeChord = " << BladeChord << endl;
    Info << "BladeTwist = " << BladeTwist << endl;
    Info << "AirfoilTypesDistinct = " << airfoilTypesDistinct << endl;
    Info << "BladeAirfoilTypeID = " << BladeAirfoilTypeID << endl << endl << endl;

    Info << "airfoilAlpha = " << airfoilAlpha << endl;
    Info << "airfoilCl = " << airfoilCl << endl;
    Info << "airfoilCd = " << airfoilCd << endl;

    Info << "influenceCells = " << influenceCells << endl << endl << endl;

    Info << "bladeDs = " << bladeDs << endl;
    Info << "towerDs = " << towerDs << endl;
    Info << "nacelleDs = " << nacelleDs << endl;
    Info << "bladePoints = " << bladePoints << endl;
    Info << "towerPoints = " << towerPoints << endl;
    Info << "nacellePoints = " << nacellePoints << endl;
    Info << "bladePointRadius = " << bladePointRadius << endl;
    Info << "towerPointHeight = " << towerPointHeight << endl;
    Info << "towerShaftIntersect = " << towerShaftIntersect << endl;
    Info << "rotorApex = " << rotorApex << endl;
    Info << "uvShaft = " << uvShaft << endl;
    Info << "uvShaftDir = " << uvShaftDir << endl;
    Info << "uvTower = " << uvTower << endl;
    Info << "deltaNacYaw = " << deltaNacYaw << endl;
    Info << "deltaAzimuth = " << deltaAzimuth << endl;

    Info << "bladePointForce = " << bladePointForce << endl;
    Info << "nacellePointForce = " << nacellePointForce << endl;
    Info << "towerPointForce = " << towerPointForce << endl;
    Info << "bladeWindVectors = " << bladeWindVectors << endl;
    Info << "nacelleWindVector = " << nacelleWindVector << endl;
    Info << "towerWindVectors = " << towerWindVectors << endl;
    Info << "bladeAlignedVectors = " << bladeAlignedVectors << endl;
}


volVectorField& horizontalAxisWindTurbinesALM_tn::force()
{
    // Return the body force field to the solver
    return bodyForce;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbineModels
} // End namespace Foam

// ************************************************************************* //

