const labelUList& own = mesh.owner();

p_L = fvc :: interpolate(p, pos, "reconstruct(p)") ;
p_R = fvc :: interpolate(p, neg, "reconstruct(p)") ;

U_L = fvc :: interpolate(U, pos, "reconstruct(U)") & ( mesh.Sf() / mesh.magSf() ) ;
U_R = fvc :: interpolate(U, neg, "reconstruct(U)") & ( mesh.Sf() / mesh.magSf() ) ;

rho_L = fvc :: interpolate(rho, pos, "reconstruct(rho)") ;
rho_R = fvc :: interpolate(rho, neg, "reconstruct(rho)") ;

surfaceVectorField Uvec_L = fvc :: interpolate(U, pos, "reconstruct(U)") ;
surfaceVectorField Uvec_R = fvc :: interpolate(U, neg, "reconstruct(U)") ;

surfaceScalarField UU_L = Uvec_L & Uvec_L ;
surfaceScalarField UU_R = Uvec_R & Uvec_R ;


c12 = fvc :: interpolate(c) ;


M_L = U_L / c12;
M_R = U_R / c12;


surfaceScalarField duc = fvc :: interpolate(ducSensor) ;

		
forAll(own,iface)

{

if(duc[iface] > ducLevelPress)
{

	scalar mp = M_L[iface] ;
	scalar mn = M_R[iface] ;
	
	scalar pl = p_L[iface] ;
	scalar pr = p_R[iface] ;
				
	scalar bl = betapcalc( mp ) ;
	scalar br = betamcalc( mn ) ;

	scalar ul = U_L[iface] ;
	scalar ur = U_R[iface] ;

	scalar uul = UU_L[iface] ;
	scalar uur = UU_R[iface] ;

	scalar Mtilde = min(1.0, Foam :: sqrt( (uul + uur) / 2 ) / c12[iface] ) ;

	scalar Chi = (1 - Mtilde) * (1 - Mtilde) ;

	
	pave[iface] += duc[iface] * ( 0.5 * ( bl - br ) * ( pl - pr ) + ( 1.0 - Chi ) * ( bl + br - 1.0 ) * 0.5 * ( pl + pr ) );

}	
	
}	



forAll( mesh.boundaryMesh() , iPatch )

{

const polyPatch& patch = mesh.boundaryMesh()[iPatch] ;

if ((patch.type()=="processor"))

{

forAll( mesh.boundaryMesh()[iPatch] , iface )

{

if(duc.boundaryField()[iPatch][iface] > ducLevelPress )
{

        scalar mp = M_L.boundaryField()[iPatch][iface] ;
        scalar mn = M_R.boundaryField()[iPatch][iface] ;

	scalar pl = p_L.boundaryField()[iPatch][iface] ;
	scalar pr = p_R.boundaryField()[iPatch][iface] ;
	
	scalar bl = betapcalc( mp ) ;
	scalar br = betamcalc( mn ) ;

	scalar ul = U_L.boundaryField()[iPatch][iface] ;
        scalar ur = U_R.boundaryField()[iPatch][iface] ;

        scalar uul = UU_L.boundaryField()[iPatch][iface] ;
        scalar uur = UU_R.boundaryField()[iPatch][iface] ;

        scalar Mtilde = min(1.0, Foam :: sqrt( (uul + uur) / 2 ) / c12.boundaryField()[iPatch][iface] ) ;

        scalar Chi = (1 - Mtilde) * (1 - Mtilde) ;


	pave.boundaryFieldRef()[iPatch][iface] += duc.boundaryField()[iPatch][iface] * ( 0.5 * ( bl - br ) * ( pl - pr ) + ( 1.0 - Chi ) * ( bl + br - 1.0 ) * 0.5 * ( pl + pr ) ) ;
	
}

}

}

}

