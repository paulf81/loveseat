        // Compute the SGS temperature diffusivity throughout the domain on all cell faces.
        {
             // first, compute turbulent Prandtl number (sensitized to stratification)
             // requiring an SGS velocity scale.  From Moeng (1984) one can equate the
             // velocity scale 0.1*sqrt(e') to the Smagorinsky velocity scale Cs^2*delta*||S_ij||.
             // Then that equality is used in Moeng's equation 18b to get the length
             // scale ls.
             
             volScalarField velScale = turbulence->nuEff()/turbulence->delta();

             // also, a measure of stability is needed.
             dimensionedScalar tiny
             (
                  "tiny",
                  dimensionSet(0,0,-2,0,0,0,0),
                  scalar(1.0E-20)
             );
             volScalarField stability = (mag(g)/TRef) * (fvc::grad(T) & nUp);

             // and, a lengthscale is needed
             volScalarField lPos = max(min(7.6*velScale*Foam::sqrt(1.0/max(tiny,stability)),turbulence->delta()),0.0075*turbulence->delta());
             volScalarField lNeg = 0.0*lPos + turbulence->delta();
             volScalarField l = (1.0-pos(stability))*lNeg + pos(stability)*lPos;
             Prt = 1.0/(1.0 + 2.0*(l/turbulence->delta()));
        }

        //   second, use the SGS viscosity and turbulent Prandtl number to compute the turbulent
        //   diffusivity of temperature.
        kappaLES = turbulence->nuEff()/Prt;


        // Update the deviatoric part of the SGS stress tensor at interior faces
        surfNormR = -( (fvc::interpolate(turbulence->nuEff()) * mesh.magSf() * fvc::snGrad(U)) +
                       (fvc::interpolate(turbulence->nuEff() * dev(fvc::grad(U)().T())) & mesh.Sf()) );

        // Update the SGS temperature flux vector at interior faces.
        surfNormQ = -( (fvc::interpolate(kappaLES) * mesh.magSf() * fvc::snGrad(T)) );

        // Update the deviatoric part of the SGS stress tensor and the SGS temperature
        // flux vector at boundary faces.  The upper boundary stress and temperature
        // flux is set to zero.  The lower boundary uses a rough model.
        #include "SGSCorrectBoundaryConditions.H"

        // Calculate the divergence of the SGS stress tensor
        divDevR = fvc::div(surfNormR);

        // Calculate the divergence of the SGS temperature flux vector.
        divQ = fvc::div(surfNormQ);

        // Calculate the deviatoric SGS stress tensor.
        devR = turbulence->devReff();

        // Calculate the SGS temperature flux vector.
        q = - kappaLES * fvc::grad(T);

        // Calculate the SGS tke.
        k = turbulence->k();
