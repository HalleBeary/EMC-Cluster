#ifndef free_flight_scatter_HPP
#define free_flight_scatter_HPP



#include "treeStructure.hpp"
#include "drift.hpp"
#include "scatter.hpp"


//* Free Flight scatter function for a single particle */

void free_flight_scatter (geometryClass DevGeometry, ContainerClass* particles, scat_paramClass *scat_par, mat_paramClass *mat_par, FSize idxPart, int tsim, int tpulse)
{   

    static double dtau, dt2, dt3, dtp, rr; // variables to handle free flight // scatter
    static double timestep = 1e-16;
  //  static double tau = 1e-16 ; // relaxation time normalisation

   // Drift(DevGeometry, particles, mat_par, idxPart, tau); 



    
    // Fetch variable to handle particle free-flight time
    dtau = particles->getPhysicalValues(0,1)[idxPart]; 


     // Initial free-flight of the carriers //
    dt2 = (dtau < timestep ? dtau : timestep); 	
    			// min 
    Drift(DevGeometry, particles, mat_par, idxPart, dt2); 	// time to advance the trajectory //

    // Free-flight and scatter part //
    if (dtau <= timestep)
      do
      {
    	  if ( scat_par->get_maxScatMech() > 0 && particles->getPhysicalValues(1,4)[idxPart] != 0)  // If there are scattering mechanisms present

            Scatter(particles, mat_par, scat_par, idxPart); // scattering

    do { rr = Randomizer(); } while (rr <= 1e-6);  
    dt3 = -log(rr) * scat_par->get_taumax(); //Generate random flight-times. Multiply random number with maximum flight time tau 
    dtp = timestep - dtau; 	// remaining time to scatter in dt-interval //
    dt2 = (dt3 <= dtp ? dt3 : dtp); // Flight time for subsequent flight equal to dt3, if remaining time allows it.

    Drift(DevGeometry, particles, mat_par, idxPart, dt2); // free flight call again. The last time, if the while condition below is not satisfied, this sets the positions/momenta for the FMM
  
    dtau += dt3; // update times. Dtau is time of that particle has flew. dt3 is time of last flight time. 
        } while (dtau < timestep); // keep while loop going until time flew of particle is equal to timestep

      dtau -= timestep;  // this 'resets' dtau. The time left after substracting timestep will be the flight time of the first inital flight in the next round (dtau)!!

    particles->getPhysicalValues(0,1)[idxPart] = dtau; // setting back dtau to particles octree
    
 
}




#endif
