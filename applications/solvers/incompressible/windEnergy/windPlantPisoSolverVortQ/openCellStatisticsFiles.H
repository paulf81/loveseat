    // Define a set of mean variable lists
    List<scalar>     TmeanLevelsList(hLevelsTotal,0.0);
    List<vector>     UmeanLevelsList(hLevelsTotal,vector::zero);
    List<symmTensor> velFluxLevelsList(hLevelsTotal,symmTensor::zero);
    List<symmTensor> velFluxFluxLevelsList(hLevelsTotal,symmTensor::zero);
    List<vector>     tempFluxLevelsList(hLevelsTotal,vector::zero);



    // Create the path to the averaging directory.
    fileName rootDir;
    if (Pstream::parRun())
    {
        rootDir = runTime.path()/"../averaging";
    }
    else
    {
        rootDir = runTime.path()/"averaging";
    }

    // If the averaging directory doesn't exist, make it.
    if (Pstream::master() && !isDir(rootDir))
    {
        mkDir(rootDir);
    }

    // Check to see if the start time directory exists within the averaging directory; if not, create it.
    if (Pstream::master() && !isDir(rootDir/runTime.timeName()))
    {
        mkDir(rootDir/runTime.timeName());
    }



    // Open and write out the cell center vertical levels file
    OFstream hLevelsFile
        (
            rootDir/runTime.timeName()/"hLevelsCell"
        );
    forAll(hLevelsValues,hLevelsI)
    {
        hLevelsFile << hLevelsValues[hLevelsI] << " ";
    }
    hLevelsFile << endl;



    // Open the statistics files.
    OFstream TmeanFile
        (
           rootDir/runTime.timeName()/"T_mean"
        );



    OFstream umeanFile
        (
            rootDir/runTime.timeName()/"U_mean"
        );
    OFstream vmeanFile
        (
            rootDir/runTime.timeName()/"V_mean"
        );
    OFstream wmeanFile
        (
            rootDir/runTime.timeName()/"W_mean"
        );



    OFstream uuPmeanFile
        (
            rootDir/runTime.timeName()/"uu_mean"
        );
    OFstream uvPmeanFile
        (
            rootDir/runTime.timeName()/"uv_mean"
        );
    OFstream uwPmeanFile
        (
            rootDir/runTime.timeName()/"uw_mean"
        );
    OFstream vvPmeanFile
        (
            rootDir/runTime.timeName()/"vv_mean"
        );
    OFstream vwPmeanFile
        (
            rootDir/runTime.timeName()/"vw_mean"
        );
    OFstream wwPmeanFile
        (
            rootDir/runTime.timeName()/"ww_mean"
        );



    OFstream wuuPmeanFile
        (
            rootDir/runTime.timeName()/"wuu_mean"
        );
    OFstream wuvPmeanFile
        (
            rootDir/runTime.timeName()/"wuv_mean"
        );
    OFstream wuwPmeanFile
        (
            rootDir/runTime.timeName()/"wuw_mean"
        );
    OFstream wvvPmeanFile
        (
            rootDir/runTime.timeName()/"wvv_mean"
        );
    OFstream wvwPmeanFile
        (
            rootDir/runTime.timeName()/"wvw_mean"
        );
    OFstream wwwPmeanFile
        (
            rootDir/runTime.timeName()/"www_mean"
        );



    OFstream TuPmeanFile
        (
            rootDir/runTime.timeName()/"Tu_mean"
        );
    OFstream TvPmeanFile
        (
            rootDir/runTime.timeName()/"Tv_mean"
        );
    OFstream TwPmeanFile
        (
            rootDir/runTime.timeName()/"Tw_mean"
        );
