    // Create and read in initial time-averaged velocity field.
    Info << "Creating time-average velocity field, UAvg..." << endl;
    volVectorField UAvg
    (
        IOobject
        (
            "UAvg",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("UAvg",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    );

    // Create and read in initial time-averaged pressure field.
    Info << "Creating time-averaged pressure field, pAvg..." << endl;
    volScalarField pAvg
    (
        IOobject
        (
            "pAvg",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pAvg",dimensionSet(0, 2, -2, 0, 0, 0, 0),0.0)
    );

    // Create and read in initial RMS velocity field.
    Info << "Creating RMS velocity field, uRMS..." << endl;
    volVectorField uRMS
    (
        IOobject
        (
            "uRMS",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("uRMS",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    );

    // Create and read in initial RMS pressure field.
    Info << "Creating RMS pressure field, pRMS..." << endl;
    volScalarField pRMS
    (
        IOobject
        (
            "pRMS",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pRMS",dimensionSet(0, 2, -2, 0, 0, 0, 0),0.0)
    );

    // Create the time-averaged Reynolds stress tensor field.
    Info<< "Creating the time-averaged Reynolds stress field, uuPrime2..." << endl;
    volSymmTensorField uuPrime2
    (
        IOobject
        (
            "uuPrime2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("uuPrime2",dimensionSet(0, 2, -2, 0, 0, 0, 0),symmTensor::zero)
    );

    // Create the time-averaged velocity-pressure correlation field.
    Info<< "Creating the time-averaged velocity-pressure correlation field, upPrime2..." << endl;
    volVectorField upPrime2
    (
        IOobject
        (
            "upPrime2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("upPrime2",dimensionSet(0, 3, -3, 0, 0, 0, 0),vector::zero)
    );

    // Create the time-averaged pressure-pressure correlation field.
    Info<< "Creating the time-averaged pressure-pressure correlation field, ppPrime2..." << endl;
    volScalarField ppPrime2
    (
        IOobject
        (
            "ppPrime2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ppPrime2",dimensionSet(0, 4, -4, 0, 0, 0, 0),0.0)
    );

    // Create the time-averaged skewness field.
    Info<< "Creating the time-averaged skewness field, uuuPrime3..." << endl;
    volSymmTensorField uuuPrime3
    (
        IOobject
        (
            "uuuPrime3",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("uuuPrime3",dimensionSet(0, 3, -3, 0, 0, 0, 0),symmTensor::zero)
    );

    // Create the time-averaged skewness field.
    Info<< "Creating the time-averaged skewness field, vuuPrime3..." << endl;
    volSymmTensorField vuuPrime3
    (
        IOobject
        (
            "vuuPrime3",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("vuuPrime3",dimensionSet(0, 3, -3, 0, 0, 0, 0),symmTensor::zero)
    );

    // Create the time-averaged skewness field.
    Info<< "Creating the time-averaged skewness field, wuuPrime3..." << endl;
    volSymmTensorField wuuPrime3
    (
        IOobject
        (
            "wuuPrime3",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("wuuPrime3",dimensionSet(0, 3, -3, 0, 0, 0, 0),symmTensor::zero)
    );


    // Declare total averaging times for average fields and correlation fields.
    scalar avgTimeSum = 0.0;
    scalar corrTimeSum = 0.0;



    // Create the path to the averaging directory.
    fileName rootDir2;
    if (Pstream::parRun())
    {
        rootDir2 = runTime.path()/"../timeAveraging";
    }
    else
    {
        rootDir2 = runTime.path()/"timeAveraging";
    }

    // If the averaging directory doesn't exist, make it.
    if (Pstream::master() && !isDir(rootDir2))
    {
        mkDir(rootDir2);
    }



    // Open and write out the cell center vertical levels file.
    OFstream hLevelsFile2
    (
        rootDir2/"hLevelsCell"
    );

    forAll(hLevelsValues,hLevelsI)
    {
        hLevelsFile2 << hLevelsValues[hLevelsI] << " ";
    }
    hLevelsFile2 << endl;


 
    // Open the statistics files.
    OFstream UMeanFile2
    (
        rootDir2/"U_mean"
    );
    OFstream VMeanFile2
    (
        rootDir2/"V_mean"
    );
    OFstream WMeanFile2
    (
        rootDir2/"W_mean"
    );
    OFstream pMeanFile2
    (
        rootDir2/"p_mean"
    );



    OFstream uRMSFile2
    (
        rootDir2/"u_RMS"
    );
    OFstream vRMSFile2
    (
        rootDir2/"v_RMS"
    );
    OFstream wRMSFile2
    (
        rootDir2/"w_RMS"
    );
    OFstream pRMSFile2
    (
        rootDir2/"p_RMS"
    );



    OFstream uuPrime2File2
    (
        rootDir2/"uu_prime2"
    );
    OFstream uvPrime2File2
    (
        rootDir2/"uv_prime2"
    );
    OFstream uwPrime2File2
    (
        rootDir2/"uw_prime2"
    );
    OFstream vvPrime2File2
    (
        rootDir2/"vv_prime2"
    );
    OFstream vwPrime2File2
    (
        rootDir2/"vw_prime2"
    );
    OFstream wwPrime2File2
    (
        rootDir2/"ww_prime2"
    );
    OFstream upPrime2File2
    (
        rootDir2/"up_prime2"
    );
    OFstream vpPrime2File2
    (
        rootDir2/"vp_prime2"
    );
    OFstream wpPrime2File2
    (
        rootDir2/"wp_prime2"
    );



    OFstream uuuPrime3File2
    (
        rootDir2/"uuu_prime3"
    );
    OFstream uuvPrime3File2
    (
        rootDir2/"uuv_prime3"
    );
    OFstream uuwPrime3File2
    (
        rootDir2/"uuw_prime3"
    );
    OFstream uvvPrime3File2
    (
        rootDir2/"uvv_prime3"
    );
    OFstream uvwPrime3File2
    (
        rootDir2/"uvw_prime3"
    );
    OFstream uwwPrime3File2
    (
        rootDir2/"uww_prime3"
    );

    OFstream vuuPrime3File2
    (
        rootDir2/"vuu_prime3"
    );
    OFstream vuvPrime3File2
    (
        rootDir2/"vuv_prime3"
    );
    OFstream vuwPrime3File2
    (
        rootDir2/"vuw_prime3"
    );
    OFstream vvvPrime3File2
    (
        rootDir2/"vvv_prime3"
    );
    OFstream vvwPrime3File2
    (
        rootDir2/"vvw_prime3"
    );
    OFstream vwwPrime3File2
    (
        rootDir2/"vww_prime3"
    );

    OFstream wuuPrime3File2
    (
        rootDir2/"wuu_prime3"
    );
    OFstream wuvPrime3File2
    (
        rootDir2/"wuv_prime3"
    );
    OFstream wuwPrime3File2
    (
        rootDir2/"wuw_prime3"
    );
    OFstream wvvPrime3File2
    (
        rootDir2/"wvv_prime3"
    );
    OFstream wvwPrime3File2
    (
        rootDir2/"wvw_prime3"
    );
    OFstream wwwPrime3File2
    (
        rootDir2/"www_prime3"
    );
