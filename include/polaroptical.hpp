#ifndef polar_HPP
#define polar_HPP

#include <stdio.h>
#include <math.h>
#include "material_paramaters.hpp"


/********************************************************************/
//  Calculation of polar optical phonon scattering rate. INelastic 
//  process >> absorption and emission. 
             
/********************************************************************/
// Polar Optical Phonon Scattering   //  page 47 book of Jacoboni.

int polar(scat_paramClass *scat_par, mat_paramClass *mat_par, int *iCount, double w0) 

{
  int i;
  double PolarConst, Frohlich, gamma, gammaprime, factor, Fzero, A, B, C, Phononnumb, ee, eeprime, polarab, polarem, polarabsorption, polaremission;
  
  FILE *output15; 
  FILE *output16; 


  //--- Setting of desired parameters //
  
//  mat_par->set_eps(17);
//  mat_par->set_eps_infty(100);
  // -------------------------------- //


  Phononnumb = 1.0 / (exp(mat_par->get_hbar() *  w0  / (mat_par->get_kb() * mat_par->get_temp())) - 1.0); 
  Frohlich   = (1.0 / mat_par->get_eps_0()) * (1.0 / mat_par->get_eps_infty() - 1.0 / mat_par->get_eps());   // in units of eps_0.
  PolarConst =  Frohlich * mat_par->get_q() * mat_par->get_q() *  sqrt(mat_par->get_effmassX() / 2.0 / mat_par->get_q()) * w0 / ( 4* M_PI * mat_par->get_hbar()) ; 


  //=== ABSORPTION ===
  
  
  polarab = PolarConst * Phononnumb; 

  output15=fopen("ratePolarAbs.csv", "w");

  ++(*iCount);
  for (i = 1; i <= NLEV; ++i) // Create scatrate for different energy levels
  {
    ee = i * scat_par->get_deltaE(); 
    eeprime = ee + mat_par->get_hbar() * w0 / mat_par->get_q();  // in eV
    gamma = ee*(1.0 + mat_par->get_alpha()*ee);             
    gammaprime = eeprime*(1.0 + mat_par->get_alpha()*eeprime);
    A = pow((2.0*(1.0 + mat_par->get_alpha()*2.0*ee)*(1.0 + mat_par->get_alpha()*eeprime) + mat_par->get_alpha()*(gamma + gammaprime)), 2);
    B = -2.0*mat_par->get_alpha()*sqrt(gamma)*sqrt(gammaprime)*(4.0*(1.0 + mat_par->get_alpha()*ee)*(1.0 + mat_par->get_alpha()*eeprime) + mat_par->get_alpha()*(gamma + gammaprime));
    C = 4.0*(1.0 + mat_par->get_alpha()*ee)*(1.0 + mat_par->get_alpha()*eeprime)*(1.0 + mat_par->get_alpha()*2.0*ee)*(1 + mat_par->get_alpha()*2.0*eeprime);
    Fzero = (1 / C) * (A*log(abs((sqrt(gamma) + sqrt(gammaprime)) / (sqrt(gamma) - sqrt(gammaprime)))) + B); 

    factor = (1.0 + 2.0 * mat_par->get_alpha() * eeprime )/ (sqrt(gamma));
//    std::cout << "gamma  " <<  gamma   << "\n\n" << "\n";
    

    polarabsorption = polarab * factor * Fzero; // quantitiy to use in scattering table

   // std::cout << "Polarabsoprtion Rate: " << polarabsorption  << "\n\n" << "\n";

    scat_par->set_scatTable(i-1, iCount, polarabsorption);

    fprintf(output15, "%f \t %e \n", ee, polarabsorption); 
  }
  fclose(output15);

  scat_par->set_flagMech(iCount, 2); // Every scattering mechanism is indexed by its counter. The Setting flag of polar phonon scattering mechanism equal to 2. 
  scat_par->set_w(iCount, w0);       // The energy won/lost for every scat. mech is given by w0. 
  scat_par->set_maxScatMech(iCount);
//  scatpar->w[*iCount - 1][*iReg - 1] =  scat_par.set  // change in energy . In scattering table use: w = scat_par.get_polarw0()

//=== EMISSION ===

  ++(*iCount); //! The Counter is added by one for emission. 

  output16=fopen("ratePolarEms.csv", "w");

  polarem = PolarConst * (Phononnumb + 1); 

  for (i = 1; i <= NLEV; ++i)
  {

    ee = i * scat_par->get_deltaE(); 
    eeprime = ee - mat_par->get_hbar() * w0 / mat_par->get_q(); // in eV -- minus sign, energy is lost


    if (eeprime <= 0)
      polaremission = 0;
      
    else
    {
      gamma = ee*(1.0 + mat_par->get_alpha()*ee);             
      gammaprime = eeprime*(1.0 + mat_par->get_alpha()*eeprime);
      A = pow((2.0*(1.0 + mat_par->get_alpha()*2.0*ee)*(1.0 + mat_par->get_alpha()*eeprime) + mat_par->get_alpha()*(gamma + gammaprime)), 2);
      B = - 2.0*mat_par->get_alpha()*sqrt(gamma)*sqrt(gammaprime)*(4.0*(1.0 + mat_par->get_alpha()*ee)*(1.0 + mat_par->get_alpha()*eeprime) + mat_par->get_alpha()*(gamma + gammaprime));
      C = 4.0*(1.0 + mat_par->get_alpha()*ee)*(1.0 + mat_par->get_alpha()*eeprime)*(1.0 + mat_par->get_alpha()*2.0*ee)*(1 + mat_par->get_alpha()*2.0*eeprime);
      Fzero = (1 / C) * (A*log(abs((sqrt(gamma) + sqrt(gammaprime)) / (sqrt(gamma) - sqrt(gammaprime)))) + B); 


      factor = (1.0 + 2.0 * mat_par->get_alpha() * eeprime )/ (sqrt(gamma));

      polaremission = polarem * factor * Fzero;
    }
    


    scat_par->set_scatTable(i-1, iCount, polaremission);

    fprintf(output16, "%f \t %e \n", ee, polaremission); 
  }
  fclose(output16);

  scat_par->set_flagMech(iCount, 2); // SAME FLAG AS ABSORPTION!! = CORRECT. +/- Energy is set by set_w function. Different iCount though! Determines different scat. probabilities of abs/ems in scat table. Here only flag is set.
  scat_par->set_w(iCount, -w0);  // - minus sign as now one phonon is emitted, energy is lowered.
  scat_par->set_maxScatMech(iCount);

std::cout << "end polar,  the total scat mechs are: " << scat_par->get_maxScatMech() << "\n\n" << "\n";
// 
  return 0;
}


#endif
