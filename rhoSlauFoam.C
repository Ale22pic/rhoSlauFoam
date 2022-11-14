//  =========                 |
//  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
//   \\    /   O peration     |
//    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
//     \\/     M anipulation  |
//-------------------------------------------------------------------------------
//License
//
//    This file is part of OpenFOAM.
//
//    OpenFOAM is free software: you can redistribute it and/or modify it
//    under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
//    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//    for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
//
//Application
//
//    rhoEnergyFoamKITAMURA 
//
//Description
//
//    Implementation of the SLAU scheme proposed by Shima and Kitamura.
//
//---------------------------------------------------------------------------

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include <fstream>      // std::ofstream

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// SLAU useful functions

   float betapcalc ( float mplus )

   {
                float betap ;

                if (abs( mplus ) < 1. )
                {
                        //betap = 0.25 * ( 2.0 - mplus ) * (mplus + 1.0 ) * (mplus + 1.0 ) ;
                        betap = 0.25 * ( 2.0 - mplus ) * Foam :: pow( (mplus + 1.0 ) , 2 ) ;
                }
                else
                {
                        betap = 0.5 * ( 1 + mplus/abs( mplus ) ) ;
                }

                return betap ;
    }

    float betamcalc ( float mminus )

    {
                float betam ;

                if (abs( mminus ) < 1. )
                {
                        //betam = 0.25 * ( 2.0 + mminus ) * (mminus - 1.0 ) * (mminus - 1.0 ) ;
                        betam = 0.25 * ( 2.0 + mminus ) * Foam :: pow( (mminus - 1.0 ) , 2 )  ;
                }
                else
                {
                        betam = 0.5 * ( 1 - mminus / abs( mminus ) ) ;
                }

                return betam ;
    }


// Main


int main(int argc, char *argv[])

// Here we start to include the files necessary for the
// solver. Descriptions for each one is given below.
 

{

    #define NO_CONTROL
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    #include "readThermophysicalProperties.H"
    #include "variables.H"


//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



    turbulence->validate() ;



//  Info message about the beginning of the loop

    
    Info << nl << "Starting time loop" << endl ;

    
//  Info message about the starting time


    Info << "Start Timing = " << runTime.clockTimeIncrement() << " s" << nl << endl ;


//  Initialize the counter for the iterations


    int Iter = 0 ;

    
//  Start time loop

    
    while (runTime.run()) 
    

    {
     

// 	Iteration counter
	

        Iter++ ; // Add 1 to the variable Iter
        

	Info << "Iter No. = " << Iter << endl ; // Print the variable Iter


//	Give info message abount the time of the current cycle using the 
//	function "timeName()" of the class "runTime"


        Info << "Iteration starting time = " << runTime.timeName() <<  endl ;
        

	Info << "Current iteration deltaT = " <<  runTime.deltaTValue() << endl;

	
	runTime++ ;


//      Saving quantities at preavious time step, those are resolved 
//      in the Runge-Kutta (RK) of 4-th order cycle.



	rhoOld = rho ; 
       
	
	rhoUOld = rhoU ; 
       
	
	rhoEOld = rhoE ; 



	for (int cycle =0 ; cycle < rkCoeff.size() ; cycle++)


	{

        

		volScalarField c = Foam :: sqrt ( 1.4 / psi ) ;



        	surfaceScalarField rhoave = fvc :: interpolate ( rho ) ;



		surfaceVectorField Uave = fvc :: interpolate ( U ) ;


 
        	surfaceScalarField phi = Uave & mesh.Sf() ;     



        	surfaceScalarField phit = fvc :: interpolate ( rhoU ) & mesh.Sf() ; 



//      	Enthalpy calculation. It is a volVectorField, so a file
//      	that containt a value of enthalpy for each cell.



		H = ( rhoE + p ) / rho ;



//      	Enthalpy at the intercell, in order to bring the
//      	enthalpy at the intercell and create a 
//      	surfaceScalarField.



        	surfaceScalarField Have = fvc :: interpolate ( H ) ;



//      	Pressure at the intercell, so from a volVectorField
//      	we pass to a surfaceScalarField, having the pressure
//      	at the intercell.



        	surfaceScalarField pave = fvc :: interpolate ( p ) ;    


		surfaceScalarField p_L = fvc :: interpolate(p, pos, "reconstruct(p)") ;       	
		surfaceScalarField p_R = fvc :: interpolate(p, neg, "reconstruct(p)") ;

       		surfaceScalarField rho_L = fvc :: interpolate(rho, pos, "reconstruct(rho)") ;
       		surfaceScalarField rho_R = fvc :: interpolate(rho, neg, "reconstruct(rho)") ;

		surfaceVectorField Uvec_L = fvc :: interpolate(U, pos, "reconstruct(U)") ;
		surfaceVectorField Uvec_R = fvc :: interpolate(U, neg, "reconstruct(U)") ;

		surfaceScalarField U_L = Uvec_L & (mesh.Sf() / (mesh.magSf())) ;
        	surfaceScalarField U_R = Uvec_R & (mesh.Sf() / (mesh.magSf())) ;

		surfaceScalarField UU_L = Uvec_L & Uvec_L ;
		surfaceScalarField UU_R = Uvec_R & Uvec_R ;
	
		surfaceScalarField H_L = fvc :: interpolate( H , pos , "reconstruct(H)" ) ;
	        surfaceScalarField H_R = fvc :: interpolate( H , neg , "reconstruct(H)" ) ;

		surfaceScalarField c12 = fvc :: interpolate(c) ;

		surfaceScalarField M_L = U_L / c12 ;
		surfaceScalarField M_R = U_R / c12 ;

        	//#include "PTerms.H"       
					      
		const labelUList& own = mesh.owner();

        	
		
		forAll(own,iface)

        	{

                scalar mp = M_L[iface] ;
                scalar mn = M_R[iface] ;

                if ( abs( mp ) < 1.0 )
                {
                        Beta_L[iface] = 0.25 * ( 2.0 - mp ) * ( mp + 1 ) * ( mp + 1 ) ;
                }
                else
                {
                        Beta_L[iface] = 0.5 * ( 1 + mp / abs( mp ) ) ;
                }

                if ( abs( mn ) < 1.0 )
                {
                        Beta_R[iface] = 0.25 * ( 2.0 + mn ) * ( mn - 1 ) * ( mn - 1 ) ;
                }
                else
                {
                        Beta_R[iface] = 0.5 * ( 1 - mn / abs( mn ) ) ;
                }
        	
		}	



		forAll( mesh.boundaryMesh() , iPatch )

        	{
                const polyPatch& patch = mesh.boundaryMesh()[iPatch] ;

                if ((patch.type()=="processor"))

                {

                forAll( patch , iface )

                {
                        scalar mp = M_L.boundaryField()[iPatch][iface] ;
                        scalar mn = M_R.boundaryField()[iPatch][iface] ;

                        if ( abs( mp ) < 1.0 )
                        {
                        Beta_L.boundaryFieldRef()[iPatch][iface] = 0.25 * ( 2.0 - mp ) * ( mp + 1 ) * ( mp + 1 ) ;
                        }
                        else
                        {
                        Beta_L.boundaryFieldRef()[iPatch][iface] = 0.5 * ( 1 + mp / abs( mp ) ) ;
                        }

                        if ( abs( mn ) < 1.0 )
                        {
                        Beta_R.boundaryFieldRef()[iPatch][iface] = 0.25 * ( 2.0 + mn ) * ( mn - 1 ) * ( mn - 1 ) ;
                        }
                        else
                        {
                        Beta_R.boundaryFieldRef()[iPatch][iface] = 0.5 * ( 1 - mn / abs( mn ) ) ;
                        }
                }

                }

        	}


		const dimensionedScalar v_zero(dimVolume/dimTime, Zero);
		const dimensionedScalar zeroVel(dimLength/dimTime, Zero);
		const dimensionedScalar zeros(dimLength/dimLength, Zero);


		surfaceScalarField Mtilde = min ( 1.0 , ( 1.0 / c12 ) * Foam :: sqrt ( 0.5 * ( UU_L ) + ( UU_R ) ) ) ;


		surfaceScalarField Chi = ( 1.0 - Mtilde ) * ( 1.0 - Mtilde ) ;


		surfaceScalarField VnMag_hat = ( rho_L * mag( U_L ) + rho_R * mag( U_R ) ) / ( rho_L + rho_R ) ;

		VnMag_hat.setOriented(false) ;

		U_L.setOriented(false) ;
		U_R.setOriented(false) ;

		surfaceScalarField fg = - max ( min ( M_L , zeros ) , neg ) * min ( max ( M_R , zeros ) , pos ) ;


		surfaceScalarField mdot = 0.5 * ( ( rho_L * U_L ) + ( rho_R * U_R ) - VnMag_hat * ( rho_R - rho_L ) ) * ( 1 - fg ) - 0.5 * Chi / c12 * ( p_R - p_L ) ;



//      	Evaluate viscous terms



		volScalarField muEff( turbulence -> muEff() ) ;



                volTensorField tauMC( "tauMC" , muEff * dev2 (Foam :: T ( fvc :: grad (U) ) ) ) ;



		surfaceScalarField muave = fvc :: interpolate ( muEff ) ; 



                volScalarField k( "k" , thermo.Cp() * muEff / Pr ) ; 



//		Then here from the center cell k is interpolated to
//		the intercell.



        	surfaceScalarField kave = fvc :: interpolate ( k ) ;



//        	Momentum viscous flux 



		surfaceVectorField momVisFlux = muave * ( fvc :: snGrad(U) * mesh.magSf() ) ; 
			       	


//		Energy viscous flux. The functions here are the same used upper,
//		here the snGrad is applied on the temperature T, which is a scalar.
  


		surfaceScalarField heatFlux =  kave * fvc :: snGrad(T) * mesh.magSf() ;



//		The momVisFlux multiplied for Uave (which is a surfaceVectorField) 
//		doing a scalar product. So the scalar product between two 
//		surfaceVectorField give a scalar, so it will give a surfaceScalarField,
//		which is the viscous work. 



	        surfaceScalarField visWork = ( momVisFlux + 
				
				fvc :: dotInterpolate(mesh.Sf(), tauMC) ) & Uave ;



//		Here is calculated the energy by summing the heat flux and
//		the viscous work.



		enVisFlux = heatFlux + visWork ;



//		Total fluxes, Eulerian - Viscous. First is calculated
//		the flux of rho by multipling the rho at the intercell
//		by phi, which is the flux at the intercell (is given
//		by the multiplication between mesh.Sf() and Uave 
//		which is the velocity vector at the interface.  



		rhoFlux = 0.5 * ( ( mdot + mag( mdot ) ) + ( mdot - mag( mdot ) ) ) * mesh.magSf() ;     	
		
		surfaceVectorField pflux = pave * mesh.Sf() ;
		//pflux.setOriented(false);

                //momVisFlux.setOriented(false);
		
		//Uvec_L.setOriented(true) ;
		//Uvec_R.setOriented(true) ;

		surfaceVectorField flmom = 0.5 * ( ( mdot + mag( mdot ) ) * Uvec_L + ( mdot - mag( mdot ) ) * Uvec_R ) * mesh.magSf() ;
	   	flmom.setOriented(true) ;	
	

		momFlux = flmom + pflux - momVisFlux ;     


		//enVisFlux.setOriented(false) ;

        	surfaceScalarField flh = 0.5 * ( ( mdot + mag( mdot ) ) * H_L + ( mdot - mag( mdot ) ) * H_R ) * mesh.magSf() ;
	       flh.setOriented(true) ;

        	enFlux = flh - enVisFlux ; 


//		Include Convective terms
		
		
		
		//#include "CTerms.H"    
						 
						
  
//		After the calculation of the flux, the divergence is
//		applied on them in order to have the elements for
//		the resolution of the RK sub-step. 



	        volScalarField rhoFl = fvc :: div ( rhoFlux ) ;

        	volVectorField momFl = fvc :: div ( momFlux ) - fvc :: div ( tauMC ) ;

                volScalarField enFl = fvc :: div ( enFlux ) ;



//      	#include "pressureGrad.H" // forcing term to keep a constant 
      	                                  // mass flow rate in channel flow.
					  // This is a file useful in 
					  // turbulent channel flow.



//         	RK sub-step are performed. 



		rho = rhoOld + rkCoeff[cycle] * runTime.deltaT() * ( -rhoFl ) ;         

 	        rhoU = rhoUOld + rkCoeff[cycle] * runTime.deltaT() * ( -momFl ) ;

	        rhoE = rhoEOld + rkCoeff[cycle] * runTime.deltaT() * ( -enFl ) ;


          
		U.ref() = rhoU() / rho() ;


//              The fields boundary condition are corrected


        	U.correctBoundaryConditions() ;

		

		rhoU.boundaryFieldRef() = rho.boundaryField() * U.boundaryField() ;

        	e = rhoE / rho - 0.5 * magSqr (U) ;
        
		e.correctBoundaryConditions();
        

        
		thermo.correct() ;


        
        	rhoE.boundaryFieldRef() = rho.boundaryField() * ( e.boundaryField() + 
				0.5 * magSqr ( U.boundaryField() ) ) ;

        
		p.ref() = rho() / psi() ;

        	p.correctBoundaryConditions() ;

	        rho.boundaryFieldRef() = psi.boundaryField() * p.boundaryField() ; // psi=1/(R*T)



        } 
        
        
	
//      End of RK time integration
	

//	"turbulence->correct()" updates the turbulence model 
//	based upon the now solved for U. What it does 
//	depends on whether you have a RANS or LES 
//	turbulence model enabled.



        turbulence -> correct() ; 
  


	runTime.write() ;
   

	
	#include "diagnostics.H" // print tke and enstrophy on diagnostics.dat



        #include "step.H"        // Evaluate Courant Number



	#include "setDeltaT.H"   // Adjust time step and it will be used at next iteration


    }

//  End of the main



    runTime.write() ;



    Info << "End\n" << endl ;



    return 0 ;

}
// ************************************************************************* //
