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

#include "wallDist.H"
#include "patchWave.H"
#include "fvMesh.H"
#include "wallPolyPatch.H"
#include "fvPatchField.H"
#include "Field.H"
#include "emptyFvPatchFields.H"
#include "wallFvPatch.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Calculate the normal vector from a three-point plane
Foam::vector Foam::wallDist::planeNormal
(
    const point& A,
    const point& B,
    const point& C
)
{
    vector normal;

    scalar vx = (B.x() - C.x());
    scalar vy = (B.y() - C.y());
    scalar vz = (B.z() - C.z());

    scalar wx = (A.x() - B.x());
    scalar wy = (A.y() - B.y());
    scalar wz = (A.z() - B.z());

    scalar vw_x = vy * wz - vz * wy;
    scalar vw_y = vz * wx - vx * wz;
    scalar vw_z = vx * wy - vy * wx;

    scalar magitude = Foam::sqrt((vw_x * vw_x) + (vw_y * vw_y) + (vw_z * vw_z));

    scalar CA = 0;
    scalar CB = 0;
    scalar CC = 0;
    if ( magitude > SMALL )
    {
        magitude = 1/magitude;

        CA = vw_x * magitude;
        CB = vw_y * magitude;
        CC = vw_z * magitude;
    }

    normal.x() = CA;
    normal.y() = CB;
    normal.z() = CC;

    return normal;
}


Foam::label Foam::wallDist::maxCompIndex(const vector& V)
{
    scalar maxComp = -VGREAT;
    label index = 0;
    forAll(V,i)
    {
        if (maxComp < mag(V[i]))
        {
            maxComp = mag(V[i]);
            index = i;
        }
    }
    return index;
}


Foam::scalar Foam::wallDist::hitDist
(
    const point& shootPt,
    const pointField& polygonPts
)
{
    vector normal = planeNormal(polygonPts[0],polygonPts[1],polygonPts[2]);

    point A = polygonPts[0];
    scalar CA = normal.x();
    scalar CB = normal.y();
    scalar CC = normal.z();
    scalar CD = 0.0 - ((CA*A.x())+(CB*A.y())+(CC*A.z()));

    scalar distance =
        mag(CA*shootPt.x()+CB*shootPt.y()+CC*shootPt.z()+CD)
      / Foam::sqrt(CA*CA+CB*CB+CC*CC);

    return distance;
}


Foam::point Foam::wallDist::hitPoint
(
    const point& shootPt,
    const pointField& polygonPts
)
{
    vector normal = planeNormal(polygonPts[0],polygonPts[1],polygonPts[2]);

    return shootPt - ((shootPt - polygonPts[0]) & normal) * normal;
}



Foam::scalar Foam::wallDist::ptInLine
(
    const point& shootPt,
    const point& es,
    const point& ee
)
{
    vector A = shootPt-es;
    vector B = ee-es;
    scalar projLen = mag(A)*(A & B)/(mag(A)*mag(B));
    if (projLen <= mag(B) && projLen > 0)
    {
        return Foam::sqrt(mag(A)*mag(A)-projLen*projLen);
    }
    else
    {
        return -1.0;
    }

}


bool Foam::wallDist::isLeft(const point& P0, const point& P1, const point& P2)
{
    if (mag(P1.x()-P0.x()) < VSMALL)
    {
        if (mag(P2.x()-P0.x()) < VSMALL)
        {
            return true;
        }

        if (P1.y() > P0.y())
        {
            if (P2.x() > P0.x())
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        else
        {
            if (P2.x() > P0.x())
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    else if (mag(P1.y()-P0.y()) < VSMALL)
    {
        if (mag(P2.y()-P0.y()) < VSMALL)
        {
            return true;
        }

        if (P1.x() > P0.x())
        {
            if (P2.y() > P0.x())
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            if (P2.y() > P0.x())
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }
    else
    {
        scalar aaa = P1.x() - P0.x();
        scalar bbb = P2.y() - P0.y();
        scalar ccc = P2.x() - P0.x();
        scalar ddd = P1.y() - P0.y();
        scalar DET = aaa*bbb-ccc*ddd;

        if (pos(DET))
        {
            return true;
        }
        else if (neg(DET))
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}


bool Foam::wallDist::wnPnPoly
(
    const point& shootPt,
    const pointField& polygonPts
)
{
    bool inPoly = 0;
    label wn = 0;
    label compDropped = maxCompIndex
    (
        planeNormal(polygonPts[0],polygonPts[1],polygonPts[2])
    );
    pointField V(polygonPts.size());
    point P;

    switch (compDropped)
    {
        case 0:
            P.x() = shootPt.y();
            P.y() = shootPt.z();
            P.z() = shootPt.x();
            break;
        case 1:
            P.x() = shootPt.x();
            P.y() = shootPt.z();
            P.z() = shootPt.y();
            break;
        case 2:
            P.x() = shootPt.x();
            P.y() = shootPt.y();
            P.z() = shootPt.z();
            break;
        default :
            WarningIn("wnPnPoly")
                << "Default is to drop the z component!"
                << endl;
            P.x() = shootPt.x();
            P.y() = shootPt.y();
            P.z() = shootPt.z();
    }

    for (label i = 0; i < V.size(); i++)
    {
        switch (compDropped)
        {
            case 0:
                V[i].x() = polygonPts[i].y();
                V[i].y() = polygonPts[i].z();
                V[i].z() = polygonPts[i].x();
                break;
            case 1:
                V[i].x() = polygonPts[i].x();
                V[i].y() = polygonPts[i].z();
                V[i].z() = polygonPts[i].y();
                break;
            case 2:
                V[i].x() = polygonPts[i].x();
                V[i].y() = polygonPts[i].y();
                V[i].z() = polygonPts[i].z();
                break;
            default :
                V[i].x() = polygonPts[i].x();
                V[i].y() = polygonPts[i].y();
                V[i].z() = polygonPts[i].z();
        }
    }

    label VN = V.size();
    V.resize(VN+1);
    V[VN] = V[0];

    for (label i = 0; i < V.size(); i++)
    {
        label pre = i;
        label nxt = i+1;
        if (nxt == V.size())
        {
            nxt = 0;
        }

        if (V[pre].y() <= P.y())
        {
            if (V[nxt].y() > P.y())
            {
                if (isLeft(V[pre], V[nxt], P))
                {
                    ++wn;
                }
            }
        }
        else
        {
            if (V[nxt].y() <= P.y())
            {
                if (!isLeft(V[pre], V[nxt], P))
                {
                    --wn;
                }
            }
        }
    }

    wn = mag(wn);
    switch (wn)
    {
        case 0: inPoly = 0; break;
        case 1: inPoly = 1; break;
        default: FatalErrorIn("wnPnPoly")
            << "How could winding number be negative? Double check please!"
            << Foam::abort(FatalError);
    }

    return inPoly;
}


Foam::label Foam::wallDist::minDistPtID
(
    const point& shootPt,
    const pointField& wallPts
)
{
    scalar minValue = VGREAT;
    label minPtID = 0;
    forAll(wallPts, pointI)
    {
        scalar dist = mag(wallPts[pointI]-shootPt);
        if (dist < minValue)
        {
            minValue = dist;
            minPtID = pointI;
        }
    }

    return minPtID;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDist::wallDist(const fvMesh& mesh, const bool correctWalls)
:
    volScalarField
    (
        IOobject
        (
            "y",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y", dimLength, GREAT)
    ),
    cellDistFuncs(mesh),
    correctWalls_(correctWalls),
    nUnset_(0)
{
    wallDist::correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDist::~wallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallDist::correct()
{
    //
    // Calculate the true distance using my own algorithm
    //
    Info<< "Using the New Wall Distance Algorithm!!" << endl;

    fvMesh mesh(cellDistFuncs::mesh());
    const fvPatchList& patches = mesh.boundary();

    List<pointField> wallPts(Pstream::nProcs());
    List<List<pointField> > wallEdges(Pstream::nProcs());
    List<List<pointField> > wallFaces(Pstream::nProcs());
    forAll(patches, patchI)
    {
        const fvPatch& currPatch = patches[patchI];

        if (isA<wallFvPatch>(currPatch))
        {
            Info<< "Wall Patch : " << currPatch.name() << endl;
            pointField wallPtsCur = currPatch.patch().localPoints();
            edgeList wallEdgesTmpList = currPatch.patch().edges();
            List<pointField> wallEdgesCur(wallEdgesTmpList.size());
            forAll(wallEdgesTmpList,edgeI)
            {
                forAll(wallEdgesTmpList[edgeI],pointI)
                {
                    wallEdgesCur[edgeI].append
                    (
                        currPatch.patch().localPoints()
                        [
                            wallEdgesTmpList[edgeI][pointI]
                        ]
                    );
                }
            }
            faceList wallFacesTmpList = currPatch.patch().localFaces();
            List<pointField> wallFacesCur(wallFacesTmpList.size());
            forAll(wallFacesTmpList,faceI)
            {
                forAll(wallFacesTmpList[faceI],pointI)
                {
                    wallFacesCur[faceI].append
                    (
                        currPatch.patch().localPoints()
                        [
                            wallFacesTmpList[faceI][pointI]
                        ]
                    );
                }
            }

            wallPts[Pstream::myProcNo()].append(wallPtsCur);
            wallEdges[Pstream::myProcNo()].append(wallEdgesCur);
            wallFaces[Pstream::myProcNo()].append(wallFacesCur);
        }
    }

    Pstream::gatherList(wallPts);
    Pstream::gatherList(wallEdges);
    Pstream::gatherList(wallFaces);
    Pstream::scatterList(wallPts);
    Pstream::scatterList(wallEdges);
    Pstream::scatterList(wallFaces);


    scalarField trueDistance(mesh.C().size(),scalar(0.0));
    forAll(mesh.C(), cellI)
    {
        point shootPt = mesh.C()[cellI];
        scalar minCCToWallDist = VGREAT;
        bool touched = 0;
        for (label procI = 0; procI < Pstream::nProcs(); ++procI)
        {
            if (wallFaces[procI].size() != 0)
            {
                forAll(wallFaces[procI],faceI)
                {
                    pointField facePts = wallFaces[procI][faceI];

                    if (wnPnPoly(hitPoint(shootPt,facePts),facePts))
                    {
                        scalar thisDist = hitDist(shootPt,facePts);
                        if (thisDist < minCCToWallDist)
                        {
                            minCCToWallDist = thisDist;
                            touched = 1;
                        }
                    }
                }
            }
        }
        if (touched)
        {
            trueDistance[cellI] = minCCToWallDist;
        }
    }

    forAll(mesh.C(), cellI)
    {
        if (trueDistance[cellI] == 0 || trueDistance[cellI] != 0)
        {
            point shootPt = mesh.C()[cellI];
            scalar minCCToWallDist = VGREAT;
            if (trueDistance[cellI] != 0)
            {
                minCCToWallDist = trueDistance[cellI];
            }
            bool touched = 0;
            for (label procI = 0; procI < Pstream::nProcs(); ++procI)
            {
                if (wallEdges[procI].size() != 0)
                {
                    forAll(wallEdges[procI],edgeI)
                    {
                        scalar tmpDist = ptInLine
                        (
                            shootPt,
                            wallEdges[procI][edgeI][0],
                            wallEdges[procI][edgeI][1]
                        );
                        if
                        (
                            (tmpDist != -1.0)
                         && (tmpDist < minCCToWallDist)
                        )
                        {
                            minCCToWallDist = tmpDist;
                            touched = 1;
                        }
                    }
                }
            }
            if (touched)
            {
                trueDistance[cellI] = minCCToWallDist;
            }
        }
    }

    forAll(mesh.C(), cellI)
    {
        if (trueDistance[cellI] == 0 || trueDistance[cellI] != 0)
        {
            point shootPt = mesh.C()[cellI];
            scalar minCCToWallDist = VGREAT;
            if (trueDistance[cellI] != 0)
            {
                minCCToWallDist = trueDistance[cellI];
            }
            for (label procI = 0; procI < Pstream::nProcs(); ++procI)
            {
                if (wallPts[procI].size() != 0)
                {
                    label minPtID = minDistPtID(shootPt,wallPts[procI]);
                    point minPt = wallPts[procI][minPtID];
                    if (mag(shootPt-minPt) < minCCToWallDist)
                    {
                        minCCToWallDist = mag(shootPt-minPt);
                    }
                }
            }
            trueDistance[cellI] = minCCToWallDist;
        }
    }

    transfer(trueDistance);


    forAll(patches, patchI)
    {
        const fvPatch& currPatch = patches[patchI];

        if (currPatch.size() != 0)
        {
            if (isA<wallFvPatch>(currPatch))
            {
                scalarField zeroField(currPatch.size(),scalar(1e-15));
                boundaryField()[patchI].transfer(zeroField);
            }
            else if (!isA<emptyFvPatch>(currPatch))
            {
                scalarField curField(currPatch.size(),scalar(0.0));

                forAll(currPatch.patch().faceCentres(), ptI)
                {
                    point shootPt = currPatch.patch().faceCentres()[ptI];

                    scalar minPatchToWallDist = VGREAT;
                    for (label procI = 0; procI < Pstream::nProcs(); ++procI)
                    {
                        if (wallPts[procI].size() != 0)
                        {
                            label minPtID = minDistPtID(shootPt,wallPts[procI]);
                            point minPt = wallPts[procI][minPtID];
                            if (mag(shootPt-minPt) < minPatchToWallDist)
                            {
                                minPatchToWallDist = mag(shootPt-minPt);
                            }
                        }
                    }

                    if (minPatchToWallDist < VSMALL)
                    {
                        curField[ptI] = -1;
                    }
                    else
                    {
                        minPatchToWallDist = VGREAT;
                        bool touched = 0;
                        for (label procI = 0; procI < Pstream::nProcs(); ++procI)
                        {
                            if (wallFaces[procI].size() != 0)
                            {
                                forAll(wallFaces[procI],faceI)
                                {
                                    pointField facePts = wallFaces[procI][faceI];

                                    if (wnPnPoly(hitPoint(shootPt,facePts),facePts))
                                    {
                                        scalar thisDist = hitDist(shootPt,facePts);
                                        if (thisDist < minPatchToWallDist)
                                        {
                                            minPatchToWallDist = thisDist;
                                            touched = 1;
                                        }
                                    }
                                }
                            }
                        }
                        if (touched)
                        {
                            curField[ptI] = minPatchToWallDist;
                        }
                    }
                }

                forAll(currPatch.patch().faceCentres(), ptI)
                {
                    if (curField[ptI] == 0 || curField[ptI] != 0)
                    {
                        point shootPt = currPatch.patch().faceCentres()[ptI];
                        scalar minPatchToWallDist = VGREAT;
                        if (curField[ptI] != 0)
                        {
                            minPatchToWallDist = curField[ptI];
                        }
                        bool touched = 0;
                        for (label procI = 0; procI < Pstream::nProcs(); ++procI)
                        {
                            if (wallEdges[procI].size() != 0)
                            {
                                forAll(wallEdges[procI],edgeI)
                                {
                                    scalar tmpDist = ptInLine
                                    (
                                        shootPt,
                                        wallEdges[procI][edgeI][0],
                                        wallEdges[procI][edgeI][1]
                                    );
                                    if
                                    (
                                        mag(tmpDist-scalar(-1.0)) > VSMALL
                                     && tmpDist < minPatchToWallDist
                                    )
                                    {
                                        minPatchToWallDist = tmpDist;
                                        touched = 1;
                                    }
                                }
                            }
                        }
                        if (touched)
                        {
                            curField[ptI] = minPatchToWallDist;
                        }
                    }
                }

                forAll(currPatch.patch().faceCentres(), ptI)
                {
                    if (curField[ptI] == 0 || curField[ptI] != 0)
                    {
                        point shootPt = currPatch.patch().faceCentres()[ptI];
                        scalar minPatchToWallDist = VGREAT;
                        if (curField[ptI] != 0)
                        {
                            minPatchToWallDist = curField[ptI];
                        }
                        for (label procI = 0; procI < Pstream::nProcs(); ++procI)
                        {
                            if (wallPts[procI].size() != 0)
                            {
                                label minPtID = minDistPtID(shootPt,wallPts[procI]);
                                point minPt = wallPts[procI][minPtID];
                                if (mag(shootPt-minPt) < minPatchToWallDist)
                                {
                                    minPatchToWallDist = mag(shootPt-minPt);
                                }
                            }
                        }
                        curField[ptI] = minPatchToWallDist;
                    }
                }

                forAll(curField,ptI)
                {
                    if (curField[ptI] == -1)
                    {
                        curField[ptI] = VSMALL;
                    }
                }

                boundaryField()[patchI].transfer(curField);
            }
        }
    }
}


// ************************************************************************* //
