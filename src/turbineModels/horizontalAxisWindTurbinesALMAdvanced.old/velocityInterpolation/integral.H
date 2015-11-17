// Use trilinear interpolation using nearest cell center and cell-center
// gradient information.
forAll(influenceCells[i], m)
{
    vector disVector = (mesh_.C()[influenceCells[i][m]] - bladePoints[i][j][k]);
    scalar dis = mag(disVector);
    if (dis <= bladeProjectionRadius[i])
    {
        // Compute the body force projection.
        scalar spreading = computeBladeProjectionFunction(disVector,i,j,k);

        // Sum up this mesh cell's contribution to the integrated velocity, 
        // weighted by the projection function.
        velocity += U_[influenceCells[i][m]] * spreading * mesh_.V()[influenceCells[i][m]];
    }
}
