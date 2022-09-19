#ifndef init_kspace_hpp
#define init_kspace_hpp

#include "treeStructure.hpp"
#include "material_paramaters.hpp"
#include "random.hpp"
#include <math.h>







void k_inreader(ContainerClass* particles, FReal* Total_Momenta, scat_paramClass *scat_par, mat_paramClass *mat_par,  FSize idxPart, int i){

    particles->getPhysicalValues(1,1)[idxPart] =  Total_Momenta[3*i] ;    
    particles->getPhysicalValues(2,1)[idxPart] =  Total_Momenta[3*i + 1] ;
    particles->getPhysicalValues(3,1)[idxPart] =  Total_Momenta[3*i + 2] ;

    if (Total_Momenta[3*i] == 0 && Total_Momenta[3*i + 1] == 0 && Total_Momenta[3*i + 2] == 0)
        {
            double k, kx, ky, kz, e, rr, fai, ct, st, tau, starting_temp;

            starting_temp = 0;

            //Get momenta 

            kx = particles->getPhysicalValues(1,1)[idxPart];    
            ky = particles->getPhysicalValues(2,1)[idxPart];
            kz = particles->getPhysicalValues(3,1)[idxPart];

            // Get Tau
            tau = particles->getPhysicalValues(0,1)[idxPart];


            // Set particle energy 
            do { rr = Randomizer(); } while (rr <= 1e-6);
            e = -  (starting_temp*1.38066e-23 / ( 1.60219e-19) * 1.5) * log(rr);// 0.5  -  (mat_par->get_vt() * 1.5) * log(rr); //energy in eV

            // Set initial free-flight //
            do { rr = Randomizer(); } while (rr <= 1e-6);
            tau = - log (rr) * scat_par->get_taumax();
        
            particles->getPhysicalValues(0,1)[idxPart] = tau;

            // Set momenta

            k = sqrt(mat_par->get_q()*2*mat_par->get_effmassE()) / (mat_par->get_hbar()) * sqrt(e * ( mat_par->get_alpha() * e + 1.0)); // factor q changes from eV to J, k in 1 / m

            fai = 2.0 * M_PI * Randomizer();


            ct  = 1.0 - Randomizer() * 2.0;
            st  = sqrt(1.0 - ct * ct);
            kx  = k * st * cos(fai);
            ky  = k * st * sin(fai);
            kz  = k * ct;

        // std::cout << "kx is: " <<  kx << "\n\n" << "\n";
        // std::cout << "ky is: " <<  ky << "\n\n" << "\n";
        // std::cout << "kz is: " <<  kz << "\n\n" << "\n";


            particles->getPhysicalValues(1,1)[idxPart] = kx;  
            particles->getPhysicalValues(2,1)[idxPart] = ky;
            particles->getPhysicalValues(3,1)[idxPart] = kz;
        }


}


void LoadParticles(EMCGenericLoader<double> *loader, const int TreeHeight, const int SubTreeHeight, OctreeClass *tree, FPoint<FReal>* Momenta, FReal* Total_Momenta){ 
  //
  FPoint<FReal> position;
  FReal physicalValue = 0.0;
  //
  FPoint<FReal>*const outParticlePositions = &position; 
  FReal*const outPhysicalValue = &physicalValue;

  for(FSize idxPart = 0 ; idxPart < loader->getNumberOfParticles() ; ++idxPart){
        FReal x,y,z;
    
    (*loader->file)  >> x >> y >> z >> (*outPhysicalValue);  
                                                            
    outParticlePositions->setPosition(x,y,z);



    if(loader->otherDataToRead> 0){
        for (FSize 	i = 0 ; i <loader->otherDataToRead; ++i){

            
            (*loader->file) >> x ;



            if (i == 2) {
              Momenta->setX(x);

              Total_Momenta[idxPart*3] = Momenta->getX();
            }
            else if (i == 3){
              Momenta->setY(x); //! Momenta Fpoint can go, is redundant

              Total_Momenta[idxPart*3 + 1] = Momenta->getY();
            }
            else if (i == 4){
              Momenta->setZ(x);

              Total_Momenta[idxPart*3 + 2] = Momenta->getZ();
            }


            
            

           
        }
    }

    tree->insert(position, idxPart, physicalValue);
  }
}

void kspace_init(ContainerClass* particles, scat_paramClass *scat_par, mat_paramClass *mat_par, FSize idxPart)
{   
    double k, kx, ky, kz, e_start, rr, fai, ct, st; // 

    //Get momenta 

    kx = particles->getPhysicalValues(1,1)[idxPart];    
    ky = particles->getPhysicalValues(2,1)[idxPart];
    kz = particles->getPhysicalValues(3,1)[idxPart];


    // Set particle energy 
    do { rr = Randomizer(); } while (rr <= 1e-6);
     e_start = - ( 300 * 1.38066e-23 ) / ( 1.60219e-19 * 1.5) * log(rr) ;

    // Set initial free-flight //
    do { rr = Randomizer(); } while (rr <= 1e-6);
     particles->getPhysicalValues(0,1)[idxPart] = - log (rr) * scat_par->get_taumax();
  
    
  //  std::cout << "energy is: " <<  e_start << "\n\n" << "\n";

    // Set momenta

    k = sqrt(mat_par->get_q()*2*mat_par->get_effmassE()) / (mat_par->get_hbar()) * sqrt(e_start * ( mat_par->get_alpha() * e_start + 1.0)); // factor q changes from eV to J, k in 1 / m

    fai = 2.0 * M_PI * Randomizer();


    ct  = 1.0 - Randomizer() * 2.0;
    st  = sqrt(1.0 - ct * ct);
    kx  = k * st * cos(fai);
    ky  = k * st * sin(fai);
    kz  = k * ct;

   // std::cout << "kx is: " <<  kx << "\n\n" << "\n";
    //std::cout << "ky is: " <<  ky << "\n\n" << "\n";
   // std::cout << "kz is: " <<  kz << "\n\n" << "\n";


    particles->getPhysicalValues(1,1)[idxPart] = kx;  
    particles->getPhysicalValues(2,1)[idxPart] = ky;
    particles->getPhysicalValues(3,1)[idxPart] = kz;

}


void deltaE(ContainerClass* particles, geometryClass DevGeometry, scat_paramClass *scat_par, mat_paramClass *mat_par, FSize idxPart)
{ 
    double pulse_energy, k, fai, ct, st, kx, ky, kz; // let the particle loop run over only fraction of total number of particles so that you only excite these particles
    
    pulse_energy = 0.3;
    
    k = sqrt(mat_par->get_q()*2*mat_par->get_effmassE()) / (mat_par->get_hbar()) * sqrt(pulse_energy * ( mat_par->get_alpha() * pulse_energy + 1.0)); // factor q changes from eV to J, k in 1 / m

    fai = 2.0 * M_PI * Randomizer();


    ct  = 1.0 - Randomizer() * 2.0;
    st  = sqrt(1.0 - ct * ct);
    kx  = k * st * cos(fai);
    ky  = k * st * sin(fai);
    kz  = k * ct;

   // std::cout << "kx is: " <<  kx << "\n\n" << "\n";
    //std::cout << "ky is: " <<  ky << "\n\n" << "\n";
   // std::cout << "kz is: " <<  kz << "\n\n" << "\n";


    particles->getPhysicalValues(1,1)[idxPart] = kx;  
    particles->getPhysicalValues(2,1)[idxPart] = ky;
    particles->getPhysicalValues(3,1)[idxPart] = kz;
    
    
}


void realspace(ContainerClass* particles, geometryClass DevGeometry, FSize idxPart, int N, int i)
{
    double posX, posY, posZ, dimfac, steps, energy;
    int l,  entries, factor;

    energy = 1.0 / 3.0 ;
    steps = pow(N, energy);
   
    dimfac = 1 / steps;

    std::vector<double> grid;

     
    for (posX = DevGeometry.Xmin + 0.5 * dimfac ; posX < DevGeometry.Xmax ; posX += dimfac ){   
        grid.push_back(posX);          
        }
   
 // -----------------------------------------------------------Grid initialization --------------------------------------------------------------

    posX = DevGeometry.Xmin + 0.5 * dimfac ;
    posY = DevGeometry.Ymin + 0.5 * dimfac ;
    posZ = DevGeometry.Zmin + 0.5 * dimfac ;

    // -------- pos X -------

    factor = round(steps);

    entries = i % factor ;
    //std::cout << " steps " << round(steps) << "\n\n" << "\n"; 
    posX = grid[entries];

    particles->getPositions()[0][idxPart] = posX;

    // -------- pos Y -------

    entries = static_cast<int>( i / steps ) ;
 
    posY = grid[entries];

    if ( static_cast<int>( i / (steps*steps)) > 0)
        {
        l = i / (steps*steps);
        posY = grid[entries - l*steps];
        }

    particles->getPositions()[1][idxPart] = posY;

    // -------- pos Z -------

    entries = static_cast<int> (i / (steps*steps)) ;


    posZ = grid[entries];

    particles->getPositions()[2][idxPart] = posZ;

   //  std::cout << "particle: " << i << ",  x " << posX << ", y: " << posY << ", z: " << posZ <<  std::endl;
// ----------------------------------------------------------- end --------------------------------------------------------------


}

void Particle_Identifier(ContainerClass* particles, FSize idxPart, int j, int N)
{
    double ratio;

    ratio = 0.5;    

    if (particles->getPhysicalValues()[idxPart] == -1)
      particles->getPhysicalValues(1,4)[idxPart] = - 1.0;   // -1 electron

    else if (particles->getPhysicalValues()[idxPart] == 1 && j <= ratio*N) // ratio determines percentage of holes vs ions and electrons
    {
        particles->getPhysicalValues(1,4)[idxPart] = 1.0; // +1 holes
    }

    else
    {
        particles->getPhysicalValues(1,4)[idxPart] = 0; // 0 ions
    }

};


void ChargeChanger(ContainerClass* particles, FSize idxPart)
{
    double dopingcharge;
    
    dopingcharge = 2.0;
    particles->getPhysicalValues()[idxPart] = dopingcharge;   

};





#endif
