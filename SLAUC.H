//if(duc.boundaryField()[iPatch][iface]> ducLevelConv)


		
forAll(own,iface)

{

if(duc[iface] > ducLevelConv)
{
	
	//rhoFlux[iface] -= mag( mdot[iface] ) * mesh.magSf()[iface] ;
	
	momFlux[iface].component(0) -= 0.5 * mag( mdot[iface] ) * ( Uvec_L[iface].component(0) - Uvec_R[iface].component(0) ) * mesh.magSf()[iface] ;
	
	momFlux[iface].component(1) -= 0.5 * mag( mdot[iface] ) * ( Uvec_L[iface].component(1) - Uvec_R[iface].component(1) ) * mesh.magSf()[iface] ;
	
	momFlux[iface].component(2) -= 0.5 * mag( mdot[iface] ) * ( Uvec_L[iface].component(2) - Uvec_R[iface].component(2) ) * mesh.magSf()[iface] ;
	
	
	enFlux[iface] -= 0.5 * mag( mdot[iface] ) * ( H_L[iface] - H_R[iface] ) * mesh.magSf()[iface] ;	
}	
	
}	



forAll( mesh.boundaryMesh() , iPatch )

{

const polyPatch& patch = mesh.boundaryMesh()[iPatch] ;

if ((patch.type()=="processor"))

{

forAll( mesh.boundaryMesh()[iPatch] , iface )

{

if(duc.boundaryField()[iPatch][iface]> ducLevelConv)

{

	//rhoFlux.boundaryFieldRef()[iPatch][iface] -= mag( mdot.boundaryField()[iPatch][iface] ) * mesh.magSf().boundaryField()[iPatch][iface] ;
	

	momFlux.boundaryFieldRef()[iPatch][iface].component(0) -= 0.5 * mag( mdot.boundaryField()[iPatch][iface] ) * ( Uvec_L.boundaryField()[iPatch][iface].component(0) - Uvec_R.boundaryField()[iPatch][iface].component(0) ) * mesh.magSf().boundaryField()[iPatch][iface] ;
	momFlux.boundaryFieldRef()[iPatch][iface].component(1) -= 0.5 * mag( mdot.boundaryField()[iPatch][iface] ) * ( Uvec_L.boundaryField()[iPatch][iface].component(1) - Uvec_R.boundaryField()[iPatch][iface].component(1) ) * mesh.magSf().boundaryField()[iPatch][iface] ;
	momFlux.boundaryFieldRef()[iPatch][iface].component(2) -= 0.5 * mag( mdot.boundaryField()[iPatch][iface] ) * ( Uvec_L.boundaryField()[iPatch][iface].component(2) - Uvec_R.boundaryField()[iPatch][iface].component(2) ) * mesh.magSf().boundaryField()[iPatch][iface] ;

	
	enFlux.boundaryFieldRef()[iPatch][iface] -= 0.5 * mag( mdot.boundaryField()[iPatch][iface] ) * ( H_L.boundaryField()[iPatch][iface] - H_R.boundaryField()[iPatch][iface] ) * mesh.magSf().boundaryField()[iPatch][iface] ;	
	
}

}

}

}

