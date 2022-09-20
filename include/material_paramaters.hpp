#ifndef material_parameters_HPP
#define material_parameters_HPP

#include "treeStructure.hpp"

#define regions 1 // number of doping regions
#define NLEV 1000 // Number of energy levels in scattering table
#define MAXSC 10 // Max number of scattering mechanisms

class geometryClass
{   

    double box_dim; 
    public: 
      
        double BoxWidth = 0.;
        FPoint <double> BoxCenter = {0., 0., 0.};
        double Xmax = 0.;
        double Xmin = 0.;
        double Ymax = 0.;
        double Ymin = 0.;
        double Zmax = 0.;
        double Zmin = 0.;

      geometryClass(FFmaGenericLoader<FReal> *loader) : box_dim{1e7}
      {
        BoxWidth=loader->getBoxWidth(); 
        BoxCenter=loader->getCenterOfBox();
        Xmax=BoxCenter.getX()+BoxWidth/2;
        Xmin=BoxCenter.getX()-BoxWidth/2;
        Ymax=BoxCenter.getY()+BoxWidth/2;
        Ymin=BoxCenter.getY()-BoxWidth/2;
        Zmax=BoxCenter.getZ()+BoxWidth/2;
        Zmin=BoxCenter.getZ()-BoxWidth/2;
      }

      double get_boxdim();



};

        double geometryClass::get_boxdim()
        {
            return box_dim;
        }


class mat_paramClass
{
        const double hbar, q, eps_0, kb, m_0, alpha, v_sound, temp, vt, density;
        double  eps, eps_infty, effmass[3];


    public:
        mat_paramClass() : hbar{1.05459e-34}, q{1.60219e-19}, eps_0{8.85419e-12}, kb{1.38066e-23},  m_0{9.11e-31},  alpha{0}, v_sound{343.0}, temp{300.0}, vt{temp*kb / ( 1.60219e-19)}, density{2000}, eps{30}, eps_infty{10}, effmass{0.2*9.11e-31, 10*0.2*9.11e-31, 0.2*9.11e-31}{}
        
        void seteffectivemasses (double effmassE, double effmassH, double effmass_ions);
        double get_effmassE ();
        double get_effmassH ();
        double get_effmass_ions ();
        double get_alpha();
        double get_hbar();
        double get_q();
        double get_kb();
        double get_eps_0(); 
        double get_eps();
	void   set_eps(double epsilon);
        double get_eps_infty();
	void   set_eps_infty(double epsilon_uv);
        double get_temp();
        double get_v_sound();
        double get_vt();
        double get_m_0();
        double get_density();


};           
        
        void mat_paramClass::seteffectivemasses (double effmassE, double effmassH, double effmass_ions)    
        {
            effmass[0] = effmassE;
            effmass[1] = effmassH;
            effmass[2] = effmass_ions;

        }

        double mat_paramClass::get_effmassE () 
        {
            return effmass[0];
        } 

        double mat_paramClass::get_effmassH ()
        {
            return effmass[1];
        } 

        double mat_paramClass::get_effmass_ions ()
        {
            return effmass[2];
        } 
	
	
        double mat_paramClass::get_alpha() 
        {
            return alpha;
        } 

        double mat_paramClass::get_hbar() 
        {
            return hbar;
        } 

        double mat_paramClass::get_kb() 
        {
            return kb;
        } 
        
        double mat_paramClass::get_q() 
        {
            return q;
        } 
        
        double mat_paramClass::get_eps_0() 
        {
            return eps_0;
        } 
        
        double mat_paramClass::get_eps() 
        {
            return eps;
        } 
	
 	void mat_paramClass::set_eps(double epsilon)
        {
            eps = epsilon;
        }
        
        double mat_paramClass::get_eps_infty() 
        {
            return eps_infty;
        } 
        
        void mat_paramClass::set_eps_infty(double epsilon_uv)
        {
            eps_infty = epsilon_uv;
        }

        double mat_paramClass::get_temp() 
        {
            return temp;
        } 
        
        double mat_paramClass::get_v_sound() 
        {
            return v_sound;
        } 
        
        double mat_paramClass::get_vt() 
        {
            return vt;
        } 
        
        double mat_paramClass::get_m_0() 
        {
            return m_0;
        } 

        double mat_paramClass::get_density()
        {
            return density;
        } 
	        


class scat_paramClass
{   
        double ScatTable[NLEV][MAXSC], def_pot, w[MAXSC], de, polarw0, taumax;
        int flagMech[MAXSC],  maxScatMech, acousticscattering, opticalscattering;


    public: 
        scat_paramClass(int ac, int op) : def_pot{-2.13 * 1.6e-19}, polarw0{5e12*2*M_PI}, acousticscattering(ac),opticalscattering(op){}

        void set_w(const int * iCount, double omega);
        double get_w(const int * iCount);
        void set_def_pot(double &d_def_pot);
        double get_def_pot();
        void set_polarw0(double *d_polarw0);
        double get_polarw0();
        double get_deltaE();
        int get_maxScatMech();
        void set_maxScatMech (const int *iCount);
        double get_scatTable (int loc, int iTop);
        void set_scatTable (int const &i, const int * iCount, double scatRate);
        void set_flagMech (const int * iCount, int input);
        int get_flagMech (int Imech);
        void set_taumax(double tau);
        double get_taumax();
        void set_acoustic();
        void set_optical();
        bool scat_acoustic();
        bool scat_optical();


};


        void scat_paramClass::set_w(const int * iCount, double omega) // Change in energy due to scattering
        {
            w[*iCount] = omega;
        }

        double scat_paramClass::get_w(const int * iCount)
        {
            return w[*iCount];
        }

        void scat_paramClass::set_def_pot(double &d_def_pot)
        {
            def_pot = d_def_pot;
        }

        double scat_paramClass::get_def_pot()    
        {
            return def_pot;
        }

        void scat_paramClass::set_polarw0(double *d_polarw0)
        {
            polarw0 = *d_polarw0;
        }

        double scat_paramClass::get_polarw0()
        {
            return polarw0;
        }

        double scat_paramClass::get_deltaE()
        {
            de = (4.0 / NLEV); 
            return de;
        }

        int scat_paramClass::get_maxScatMech ()
        {
            return maxScatMech;
        }

        void scat_paramClass::set_maxScatMech (const int *iCount)
        {
            maxScatMech = *iCount;
        }
    
        double scat_paramClass::get_scatTable (int loc, int iTop) // returns the energy at the index 'loc'
        {
            return ScatTable[loc][iTop]; 
        }

        void scat_paramClass::set_scatTable (int const &i, const int *iCount, double scatRate)
        {
            ScatTable[i][*iCount-1] = scatRate;
        }
        
        void scat_paramClass::set_flagMech (const int * iCount, int input)
        {
            flagMech[*iCount-1] = input;
        }

        int scat_paramClass::get_flagMech (int Imech)
        {
            return flagMech[Imech];
        }

        void scat_paramClass::set_taumax (double tau)
        {
            taumax = tau;
        }

        double scat_paramClass::get_taumax()
        {
            return taumax;
        }
        
        void scat_paramClass::set_acoustic ()
        {
            acousticscattering = 1;
        } 

        void scat_paramClass::set_optical ()
        {
            opticalscattering = 1;
        } 

        bool scat_paramClass::scat_acoustic()
        {
          return acousticscattering == 1;   
        }

        bool scat_paramClass::scat_optical()
        {
            return opticalscattering == 1;
        }

        



#endif
