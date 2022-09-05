#ifndef scat_table_HPP
#define scat_table_HPP

#include "treeStructure.hpp"
#include "acoustic.hpp"
#include "polaroptical.hpp"
#include "material_paramaters.hpp"


// Scattering table. Scattering mechanisms included: 
//- Isotropic scattering (Acoustic phonons) -> flag_mech = 1 ( Intervalley scattering is also isotropic, so if this is added it gets the same flag_mech as it scatters in the same way.)
//- Polar Optical Phonons -> flag_mech =2
//


void RenormalizeTable(scat_paramClass *scat_par, int *iCount) 
{
    static int i, k;
    static double x, y, z;
    static int imax;
    static double tau;
//    FSize idxPart;

    scat_par->set_maxScatMech(iCount);
    imax = scat_par->get_maxScatMech();
    std::cout << "Number of scatterin mechanisms: " << imax << std::endl;


    if (imax == 0) 
    {
      scat_par->set_taumax(2e-15);
    }
    if (imax >= 1)
    {
        if (imax > 1) // setting of scattering rates for all energies and all scattering mechanisms by adding them up (cumulative)
            for (i=2 ; i <= imax; ++i) // loop over all scattering mechanisms
                for (k = 1 ; k <= NLEV; ++k) // loop over all energy lvls
                {
                     x = scat_par->get_scatTable(k-1, i-2); 
                   // std::cout << "scat rates x: " << x  << "\n\n" << "\n";
                     y = scat_par->get_scatTable(k-1, i-1);
                   // std::cout << "scat rates third column: " <<  scat_par->get_scatTable(k-1, 0)  << "\n\n" << "\n";
                     y += x;
                     
                      
                     scat_par->set_scatTable(k-1, &i, y);   // Scatering table is constructed for every scat. mech (iCount (i)) for every energy lvl ( k)
                    
                }
        tau = 0.0;  
        for (i = 1; i <= NLEV; ++i)
            if (scat_par->get_scatTable(i-1, imax -1) > tau)
                tau = scat_par->get_scatTable(i-1, imax -1);  // setting of total scattering rate, tau, for each energy lvl. Tau not dependent on i, why the for-loop?? >> !To see when (at which energy) tau is max!

        
       // std::cout << "tau: " << tau << "\n";
        for (i = 1; i <= imax; ++i)
            for (k = 1; k <= NLEV; ++k)
            {   
   
                z = scat_par->get_scatTable(k - 1, i - 1);
               // std::cout << "k: " << k  << "\n\n" << "\n";
               // std::cout << "scat rates before renormalisation: " << z  << "\n\n" << "\n";
                z /= tau;                                   // renormalising of scattering table. Each entry is divided by the total scattering rate
                scat_par->set_scatTable(k-1, &i, z);   
               // std::cout << "scat rates after renormalisation: " << z  << "\n\n" << "\n";           
            }

     

    // Set maximum scattering tau
    
    scat_par->set_taumax(1 / tau);
    }

    std::cout << "taumax: " << scat_par->get_taumax()  << "\n\n" << "\n";


    
//    for (FSize idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart)
//    {
//      particles->getPhysicalValues(0,1)[idxPart] = tau;
//    }

//    particles->getPhysicalValues(0,1)[idxPart] = tau; // setting tau to max. scattering rate.


    std::cout << "That worked out" << std::endl;
   
}

void ScatteringTable(scat_paramClass *scat_par, mat_paramClass *mat_par)
{
    static double w0; // w0 is a scattering quantity, declared here.
    static int iCount;

//    std::cout << "scatTable\n";
    iCount = 0; // The Counter is incremented at each different scattering mechanism. The flag is equal set to the iCount.

    //scat_par->set_acoustic();
   // scat_par->set_optical();

    
    if (scat_par->scat_acoustic() == 1)
    {
        std::cout << " here ? " << std::endl;
        acoustic(scat_par, mat_par, &iCount);
    }


    if (scat_par->scat_optical() == 1)
    {
       // std::cout << " or was it here ? " << std::endl;
        w0 = scat_par->get_polarw0();  // char. phonon freq.
        polar(scat_par, mat_par, &iCount, w0);
    }
    
   std::cout << "iCount = " << iCount << ", " << &iCount << "\n";
    
    

    RenormalizeTable(scat_par, &iCount); // iCount is incremented for each scattering mechanism. The total number is passed on to RenormalizeTable.
}

#endif
