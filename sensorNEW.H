// Evaluate Ducors sensor
   volScalarField divU  = fvc::div(U) ;
   volScalarField divU2 = divU*divU ;
   volScalarField rotU2 = magSqr(fvc::curl(U)) ;
   const volScalarField nu_t = turbulence->nut();
   const volScalarField k_visc = turbulence->nu();


// surfaceScalarField dx
//   (
//    1./mesh.surfaceInterpolation::deltaCoeffs()
//   );

// volVectorField UFil   = U   + fvc::laplacian(pow(dx,2)/4.,U);
// volScalarField rhoFil = rho + fvc::laplacian(pow(dx,2)/4.,rho);
   volVectorField UFil = fvc::surfaceSum
     (
     mesh.magSf()*fvc::interpolate(U)
     )/fvc::surfaceSum(mesh.magSf());
//
   volScalarField rhoFil = fvc::surfaceSum
     (
     mesh.magSf()*fvc::interpolate(rho)
     )/fvc::surfaceSum(mesh.magSf());

// Internal field
   forAll(U,icell)
   {
    scalar filtSensorRho = Foam::pow((mag(rho[icell]-rhoFil[icell]))/(rho[icell]+rhoFil[icell] + nu_t[icell]/k_visc[icell]+1e-9),0.25);
    scalar filtSensorU = (mag(U[icell].component(0)-UFil[icell].component(0)))/(mag(U[icell].component(0))+mag(UFil[icell].component(0)) + nu_t[icell]/k_visc[icell]+1e-9);
    scalar filtSensorV = (mag(U[icell].component(1)-UFil[icell].component(1)))/(mag(U[icell].component(1))+mag(UFil[icell].component(1)) + nu_t[icell]/k_visc[icell]+1e-9);
    scalar filtSensorW = (mag(U[icell].component(2)-UFil[icell].component(2)))/(mag(U[icell].component(2))+mag(UFil[icell].component(2)) + nu_t[icell]/k_visc[icell]+1e-9);

    scalar ducSens = max(-divU[icell]/Foam::sqrt(divU2[icell] + rotU2[icell] + eps[icell]),0.) ; // detect compressions

    scalar  sensor = max(filtSensorU,filtSensorV) ;
    sensor = max(sensor,filtSensorW);
    sensor = max(sensor,filtSensorRho);
    sensor = max(sensor,ducSens);
    ducSensor[icell] = sensor;
   }
//
// Boundary field
   forAll( mesh.boundaryMesh(), iPatch )
   {
    const polyPatch& patch = mesh.boundaryMesh()[iPatch] ;
    if ((patch.type()=="processor") )
     {
      forAll( patch , iface )
      {
      /*Edited*/
       scalar filtSensorRho = Foam::pow((mag(rho.boundaryField()[iPatch][iface]-rhoFil.boundaryField()[iPatch][iface]))/(mag(rho.boundaryField()[iPatch][iface])+mag(rhoFil.boundaryField()[iPatch][iface]) + nu_t.boundaryField()[iPatch][iface]/k_visc.boundaryField()[iPatch][iface]+1e-9),0.25);
       scalar filtSensorU = (mag(U.boundaryField()[iPatch][iface].component(0)-UFil.boundaryField()[iPatch][iface].component(0)))/(mag(U.boundaryField()[iPatch][iface].component(0))+mag(UFil.boundaryField()[iPatch][iface].component(0)) + nu_t.boundaryField()[iPatch][iface]/k_visc.boundaryField()[iPatch][iface]+1e-9);
       scalar filtSensorV = (mag(U.boundaryField()[iPatch][iface].component(1)-UFil.boundaryField()[iPatch][iface].component(1)))/(mag(U.boundaryField()[iPatch][iface].component(1))+mag(UFil.boundaryField()[iPatch][iface].component(1)) + nu_t.boundaryField()[iPatch][iface]/k_visc.boundaryField()[iPatch][iface]+1e-9);
       scalar filtSensorW = (mag(U.boundaryField()[iPatch][iface].component(2)-UFil.boundaryField()[iPatch][iface].component(2)))/(mag(U.boundaryField()[iPatch][iface].component(2))+mag(UFil.boundaryField()[iPatch][iface].component(2)) + nu_t.boundaryField()[iPatch][iface]/k_visc.boundaryField()[iPatch][iface]+1e-9);

       scalar ducSens = max(-divU.boundaryField()[iPatch][iface]/Foam::sqrt(divU2.boundaryField()[iPatch][iface] + rotU2.boundaryField()[iPatch][iface] + eps.boundaryField()[iPatch][iface]),0.) ; // detect compressions

       scalar  sensor = max(filtSensorU,filtSensorV) ;
       sensor = max(sensor,filtSensorW) ;
       sensor = max(sensor,filtSensorRho);
       sensor = max(sensor,ducSens) ;

       ducSensor.boundaryFieldRef()[iPatch][iface] = sensor ; //detect compressions
 // ducSensor.boundaryFieldRef()[iPatch][iface] = Foam::sqrt(divU2.boundaryField()[iPatch][iface]/(divU2.boundaryField()[iPatch][iface] + rotU2.boundaryField()[iPatch][iface] +  nu_t.boundaryField()[iPatch][iface]/k_visc.boundaryField()[iPatch][iface]+1e-6)) ; //detect compressions & expansions

      }
     }
   }

