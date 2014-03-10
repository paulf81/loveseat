    // Create the path to the averaging directory.
    fileName rootDir;
    if (Pstream::parRun())
    {
        rootDir = runTime.path()/"../horizontalAveraging";
    }
    else
    {
        rootDir = runTime.path()/"horizontalAveraging";
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



    // Open and write out the cell center vertical levels file.
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
    OFstream UMeanFile
    (
        rootDir/runTime.timeName()/"U_mean"
    );
    OFstream VMeanFile
    (
        rootDir/runTime.timeName()/"V_mean"
    );
    OFstream WMeanFile
    (
        rootDir/runTime.timeName()/"W_mean"
    );
    OFstream pMeanFile
    (
        rootDir/runTime.timeName()/"p_mean"
    );


    OFstream uRMSFile
    (
        rootDir/runTime.timeName()/"u_RMS"
    );
    OFstream vRMSFile
    (
        rootDir/runTime.timeName()/"v_RMS"
    );
    OFstream wRMSFile
    (
        rootDir/runTime.timeName()/"w_RMS"
    );
    OFstream pRMSFile
    (
        rootDir/runTime.timeName()/"p_RMS"
    );



    OFstream uuPrime2File
    (
        rootDir/runTime.timeName()/"uu_prime2"
    );
    OFstream uvPrime2File
    (
        rootDir/runTime.timeName()/"uv_prime2"
    );
    OFstream uwPrime2File
    (
        rootDir/runTime.timeName()/"uw_prime2"
    );
    OFstream vvPrime2File
    (
        rootDir/runTime.timeName()/"vv_prime2"
    );
    OFstream vwPrime2File
    (
        rootDir/runTime.timeName()/"vw_prime2"
    );
    OFstream wwPrime2File
    (
        rootDir/runTime.timeName()/"ww_prime2"
    );
    OFstream upPrime2File
    (
        rootDir/runTime.timeName()/"up_prime2"
    );
    OFstream vpPrime2File
    (
        rootDir/runTime.timeName()/"vp_prime2"
    );
    OFstream wpPrime2File
    (
        rootDir/runTime.timeName()/"wp_prime2"
    );



    OFstream uuuPrime3File
    (
        rootDir/runTime.timeName()/"uuu_prime3"
    );
    OFstream uuvPrime3File
    (
        rootDir/runTime.timeName()/"uuv_prime3"
    );
    OFstream uuwPrime3File
    (
        rootDir/runTime.timeName()/"uuw_prime3"
    );
    OFstream uvvPrime3File
    (
        rootDir/runTime.timeName()/"uvv_prime3"
    );
    OFstream uvwPrime3File
    (
        rootDir/runTime.timeName()/"uvw_prime3"
    );
    OFstream uwwPrime3File
    (
        rootDir/runTime.timeName()/"uww_prime3"
    );

    OFstream vuuPrime3File
    (
        rootDir/runTime.timeName()/"vuu_prime3"
    );
    OFstream vuvPrime3File
    (
        rootDir/runTime.timeName()/"vuv_prime3"
    );
    OFstream vuwPrime3File
    (
        rootDir/runTime.timeName()/"vuw_prime3"
    );
    OFstream vvvPrime3File
    (
        rootDir/runTime.timeName()/"vvv_prime3"
    );
    OFstream vvwPrime3File
    (
        rootDir/runTime.timeName()/"vvw_prime3"
    );
    OFstream vwwPrime3File
    (
        rootDir/runTime.timeName()/"vww_prime3"
    );

    OFstream wuuPrime3File
    (
        rootDir/runTime.timeName()/"wuu_prime3"
    );
    OFstream wuvPrime3File
    (
        rootDir/runTime.timeName()/"wuv_prime3"
    );
    OFstream wuwPrime3File
    (
        rootDir/runTime.timeName()/"wuw_prime3"
    );
    OFstream wvvPrime3File
    (
        rootDir/runTime.timeName()/"wvv_prime3"
    );
    OFstream wvwPrime3File
    (
        rootDir/runTime.timeName()/"wvw_prime3"
    );
    OFstream wwwPrime3File
    (
        rootDir/runTime.timeName()/"www_prime3"
    );
