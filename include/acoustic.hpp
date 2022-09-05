#ifndef acoustic_HPP
#define acoustic_HPP

#include "treeStructure.hpp"
#include "material_paramaters.hpp"



// this is elastic acoustic phonon interaction, both emission and absorption. (the interaction probability with acoustic phonons)
void acoustic(scat_paramClass *scat_par, mat_paramClass *mat_par, int *iCount) 

{
    int i;
    double scatRate, Const, ee, fe, cL;
    FILE *output3;
    
    // Calculate scattering constants:
    
    cL = mat_par->get_density() * mat_par->get_v_sound() * mat_par->get_v_sound();
    

    Const = sqrt(2.0 * mat_par->get_effmassX()) * mat_par->get_effmassX() * mat_par->get_kb() * mat_par->get_temp() *
            scat_par->get_def_pot() * scat_par->get_def_pot() / (M_PI * cL * pow(mat_par->get_hbar(), 4));
    
    output3 = fopen("rateAcoustic.csv", "w");
   
    //Create scattering table
    ++(*iCount); // number of active scattering mechanisms is incremented, when scat. table is constructed.

    std::cout << "delta energy: " << scat_par->get_deltaE() << std::endl;


    for (i = 1; i <= 1000; ++i)
    {
        ee = i * scat_par->get_deltaE();

        fe = ee * (mat_par->get_alpha() * ee + 1.0);

      //  std::cout << " scatRate: " << Const << std::endl;
        scatRate = Const * sqrt(fe) * (2 * mat_par->get_alpha() * ee + 1.0);

        scat_par->set_scatTable(i-1, iCount, scatRate); // setting for each energy (i), for this ( "icount'th") scat. mechanism, the scattering rate (scatrate) equal to the scattering rate (scatRate)
        fprintf(output3, "%f \t %e \n", ee, scatRate);  
    }
    fclose(output3);
     
   scat_par->set_flagMech(iCount, 1); // Setting flag of acoustic phonon scattering mechanism equal to 1
   scat_par->set_maxScatMech(iCount);
    
    
  std::cout << "end acoustic, counter is: " << *iCount << "\n\n " <<  std::endl;

}




#endif
