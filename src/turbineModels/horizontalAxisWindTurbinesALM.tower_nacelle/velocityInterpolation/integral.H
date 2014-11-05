// Use trilinear interpolation using nearest cell center and cell-center
// gradient information.
forAll(influenceCells[i], m)
{
    scalar dis = mag(mesh_.C()[influenceCells[i][m]] - point);
    if (dis <= bladeProjectionRadius[i])
    {
        // Currently only uniform Gaussian blade distribution is implemented.
        scalar spreading = 1.0;
        if (bladeForceProjectionType[i] == "uniformGaussian")
        {
            spreading = uniformGaussian(bladeEpsilon[i][0], dis);
        }
        else
        {
            spreading = uniformGaussian(bladeEpsilon[i][0], dis);
        }

        velocity += U_[influenceCells[i][m]] * spreading * mesh_.V()[influenceCells[i][m]];
    }
}
