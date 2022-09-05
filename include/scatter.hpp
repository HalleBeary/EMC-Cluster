#ifndef scatter_HPP
#define scatter_HPP

#include "treeStructure.hpp"
#include "random.hpp"
#include "math.h"

//* Different scattering mechanisms scatter in different ways:
//* Isotropic: - Acoustic phonons
//* Anisotropic: - Polar Optical Phonons

// acoustical phonon scattering is isotropic

void Isotropic (ContainerClass* particles, mat_paramClass *mat_par, FSize idxPart)
{
    static double energy, kx, ky, kz, rknew, fi, ct, st, gamma;

    kx = particles->getPhysicalValues(1,1)[idxPart];
    ky = particles->getPhysicalValues(2,1)[idxPart];
    kz = particles->getPhysicalValues(3,1)[idxPart];

    gamma   = 1 / ( 2 * mat_par->get_q() * mat_par->get_effmassX()) * mat_par->get_hbar()* mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz) ; // parabolic energy, alfa -> 0

    if (mat_par->get_alpha() > 0)
        energy  =  1 / (mat_par->get_q()) * 0.5*(1 / mat_par->get_alpha())*(pow(1 + 4*gamma*mat_par->get_alpha(), 0.5) - 1); 
    else
        energy = gamma ; //   energy = 1 / (mat_par->get_q()) * mat_par->get_hbar()*mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz)/(2*0.15*mat_par->get_m_0());
    
    rknew = /* some constants */  pow(energy*(mat_par->get_alpha() * energy + 1.0), 0.5);

    fi = 2.0 * M_PI * Randomizer();
    ct = 1.0 - Randomizer() * 2.0;
    st = pow(1.0 - ct * ct, 0.5);

    particles->getPhysicalValues(1,1)[idxPart] = rknew * st * cos(fi);
    particles->getPhysicalValues(2,1)[idxPart] = rknew * st * sin(fi);
    particles->getPhysicalValues(3,1)[idxPart] = rknew * ct;

}

// Change in momentum due to polar optical phonon scattering

void Anisotropic (ContainerClass* particles, mat_paramClass *mat_par, scat_paramClass *scat_par, FSize idxPart, int iFix)
{
    static double rr, e, eprime, ge, gnew, kx, ky, kz, kxy, k, kp, cth0, sth0, cfi0, sfi0, f, cth, sth, fi, kxp, kyp, kzp, gamma;

    kx = particles->getPhysicalValues(1,1)[idxPart];
    ky = particles->getPhysicalValues(2,1)[idxPart];
    kz = particles->getPhysicalValues(3,1)[idxPart];

    gamma   = 1 / ( 2 * mat_par->get_q() * mat_par->get_effmassX()) * mat_par->get_hbar()* mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz) ; // parabolic energy, alfa -> 0

    if (mat_par->get_alpha() > 0)
        e  =  1 / (mat_par->get_q()) * 0.5*(1 / mat_par->get_alpha())*(pow(1 + 4*gamma*mat_par->get_alpha(), 0.5) - 1); 
    else
        e = gamma ;
    
    // e = 1 / (mat_par->get_q()) * mat_par->get_hbar()*mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz)/(2.0*mat_par->get_effmassX());
    
    eprime = e + scat_par->get_w(&iFix) * mat_par->get_hbar() / mat_par->get_q();


//    std::cout <<"e:   " << e << "\n\n" << "\n"; 
 //   std::cout <<"eprime:   " << eprime << "\n\n" << "\n"; 

    if (eprime < 0) {
        eprime = e;

    }
    else {
    // Rotation Angles of k (Vasileska/Goodnick) 

        kxy  = sqrt(kx * kx + ky * ky);      
        k    = sqrt(kxy * kxy + kz * kz);      
        cth0 = kz / k; 
        sth0 = kxy / k; 
        cfi0 = kx / kxy; 
        sfi0 = ky / kxy;

        //    std::cout << " here ? " << std::endl;
        
        kp   =( 1 / mat_par->get_hbar() ) * sqrt(2.0 * mat_par->get_effmassX() * mat_par->get_q() * eprime * ( 1.0 + mat_par->get_alpha() * eprime)); // factor q for eV conversion
        
     //   std::cout << "kp: " <<  kp << "\n\n" << "\n";  

        ge   = e * (1 + mat_par->get_alpha() * e);
        gnew = eprime * (1 + mat_par->get_alpha() * eprime);
        f    = 2 * sqrt(ge * gnew) / (ge + gnew - 2*sqrt(ge * gnew));
        rr   = Randomizer();
        cth  = ((f + 1) - pow(1 + 2 * f, rr))/ f;
        sth  = sqrt(1 - cth*cth);
        fi   = 2.0 * M_PI * Randomizer();

        kxp  = kp * sth * cos(fi);
        kyp  = kp * sth * sin(fi);
        kzp  = kp * cth;

        kx   = kxp * cfi0 * cth0 - kyp * sfi0 + kzp * cfi0 * sth0;      
        ky   = kxp * sfi0 * cth0 + kyp * cfi0 +  kzp * sfi0 * sth0;      
        kz   = - kxp * sth0 + kzp * cth0;


        particles->getPhysicalValues(1,1)[idxPart] = kx;
        particles->getPhysicalValues(2,1)[idxPart] = ky;
        particles->getPhysicalValues(3,1)[idxPart] = kz;

    }
 //   std::cout <<"We went through optical phonon scattering here" ;
}

void Anisotropic_Isotropic(ContainerClass* particles, mat_paramClass *mat_par, scat_paramClass *scat_par, FSize idxPart, int iFix) // change energy, but scatter isotropically
{
    static double e, eprime, k,  kx, ky, kz, fi, ct, st, gamma;

    
    kx = particles->getPhysicalValues(1,1)[idxPart];
    ky = particles->getPhysicalValues(2,1)[idxPart];
    kz = particles->getPhysicalValues(3,1)[idxPart];


    
    gamma   = 1 / ( 2 * mat_par->get_q() * mat_par->get_effmassX()) * mat_par->get_hbar()* mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz) ; // parabolic energy, alfa -> 0

    if (mat_par->get_alpha() > 0)
        e  =  1 / (mat_par->get_q()) * 0.5*(1 / mat_par->get_alpha())*(pow(1 + 4*gamma*mat_par->get_alpha(), 0.5) - 1); 
    else
        e = gamma ;  //  e = 1 / (mat_par->get_q()) * mat_par->get_hbar()*mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz)/(2.0*mat_par->get_effmassX());


    eprime = e + scat_par->get_w(&iFix) * mat_par->get_hbar() / mat_par->get_q();


  
        
    
    if (eprime < 0) {
        eprime = e;

    }
    else {
        
        k = /* some constants */sqrt(mat_par->get_q()*2*mat_par->get_effmassX()) / (mat_par->get_hbar()) *  pow(eprime*(mat_par->get_alpha() * eprime + 1.0), 0.5);
    
        fi = 2.0 * M_PI * Randomizer();

        ct = 1.0 - Randomizer() * 2.0;
        st = pow(1.0 - ct * ct, 0.5);

        kx = k * st * cos(fi);
        ky = k * st * sin(fi);
        kz = k * ct;


        particles->getPhysicalValues(1,1)[idxPart] = kx;
        particles->getPhysicalValues(2,1)[idxPart] = ky;
        particles->getPhysicalValues(3,1)[idxPart] = kz;


    }

    
   // std::cout <<"We went through optical phonon scattering here" ;
}

void Scatter(ContainerClass* particles, mat_paramClass *mat_par, scat_paramClass *scat_par, FSize idxPart)
{
    static int i, iMech, loc, selectMech;
    static int iTop;
    static double kx, ky, kz, energy, rr, boundLower, boundUpper, gamma;   

    //get particle energy

    kx = particles->getPhysicalValues(1,1)[idxPart];
    ky = particles->getPhysicalValues(2,1)[idxPart];
    kz = particles->getPhysicalValues(3,1)[idxPart];

    gamma   = 1 / ( 2 * mat_par->get_q() * mat_par->get_effmassX()) * mat_par->get_hbar()* mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz) ; // parabolic energy, alfa -> 0

    if (mat_par->get_alpha() > 0)
        energy  =  1 / (mat_par->get_q()) * 0.5*(1 / mat_par->get_alpha())*(pow(1 + 4*gamma*mat_par->get_alpha(), 0.5) - 1); 
    else
        energy = gamma ; // energy = 1 / (mat_par->get_q()) * mat_par->get_hbar()* mat_par->get_hbar()*(kx*kx + ky*ky + kz*kz)/(2*mat_par->get_effmassX());


//     /* some constants */  kx*kx + ky*ky + kz*kz;
    
    // Calculate index to the scattering table. scat. table is 2d array. First index (energy) is determined. Then, the random number can be compared with different 
    // scattering probabilities at that energy ( they are cumulative)) 
    //std::cout <<"1" << energy << std::endl;
    

    loc = static_cast<int> (energy / scat_par->get_deltaE());

   // std::cout <<"loc" << loc << std::endl;
   // loc = (int) (energy / scat_par.get_deltaE()); // energy particle divided by delta energy per lvl, yields index to scattering table lvl.

    if (loc == 0)                 // boundary conditions
        loc = 1;
    if (loc > NLEV)
        loc = NLEV;

    // selection of scattering mechanism 

    iTop = scat_par->get_maxScatMech(); // Max. number of scattering mechanisms. This is set in scattering table.

    rr = Randomizer();
  

    if (rr >= scat_par->get_scatTable(loc - 1, iTop - 1))  // if random number is bigger then maximum probability of scattering
        return; // self-scattering
    
    if (rr < scat_par->get_scatTable(loc - 1, 0)) // # of scattering mechanisms is only 1. (iTop = 1) 

        iMech = 1;
    
    else if (iTop > 1)
        for (i = 1; i <= iTop -1; ++i) 
        {
            boundLower = scat_par->get_scatTable(loc, i - 1); // At energy index 'loc', for each scattering mechanism the lower and upper bounds are determined. 
            boundUpper = scat_par->get_scatTable(loc, i);
            if (rr >= boundLower && rr < boundUpper) // If random number is between bounds, the mechanism is selected. Imech is set to i+1, and the loop is ended.
            {
                iMech = i + 1; 
                break;
            }
        }
    

    selectMech = scat_par->get_flagMech(iMech -1); // Flags are set. Same for abs/ems polar

    //  std::cout <<"selectMech is: "<< selectMech << std::endl;
   // std::cout <<"2" << energy << std::endl;

    if      (selectMech == 1) Isotropic(particles, mat_par, idxPart); // scattering mechanism 1, isotropic etc;
  //  if      (selectMech == 1) Anisotropic(particles, mat_par, scat_par, idxPart, iMech);
    else if (selectMech == 2) Anisotropic(particles, mat_par, scat_par, idxPart, iMech); // iMech defines what type of scattering mechanism (absorption or emission)
  //  else if (selectMech == 3); // etc

    // updating momenta is done in different types of scattering functions (isotropic, polar)


}



#endif
