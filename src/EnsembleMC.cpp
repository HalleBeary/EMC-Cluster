#include "treeStructure.hpp"

#include "resetParticleLocation.hpp"
#include "EMCGenericWriter.hpp"
#include "free_flight_scatter.hpp"
#include "scat_table.hpp"
#include "acoustic.hpp"
#include "polaroptical.hpp"
#include "output.hpp"
#include "init_kspace.hpp"

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{ 
  // Define the mover and arranger classes dealing with moving and rearranging particles
  typedef FBasicParticleContainerIndexedMover<FReal,OctreeClass, ContainerClass> MoverClass;
  typedef FOctreeArranger<FReal, OctreeClass, ContainerClass, MoverClass> ArrangerClass;
  typedef FArrangerPeriodic<FReal, OctreeClass, ContainerClass, MoverClass> ArrangerClassPeriodic;

  const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZQ100.bfma" );
  const std::string filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
  const int TreeHeight       = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
  const int SubTreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
  bool periodicCondition     = false ;     
  if(FParameters::existParameter(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options))
      periodicCondition = true;
  const unsigned int aboveTree = FParameters::getValue(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options, -1);

#ifdef _OPENMP
  const unsigned int NbThreads = periodicCondition ? omp_get_max_threads() : // Setting of number of threads when algo is running periodically
                                 FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, omp_get_max_threads());
  omp_set_num_threads(NbThreads);
  std::cout << "\n>> Using " << NbThreads << " threads.\n" << std::endl;
#else
  const int NbThreads =  1;
  std::cout << "\n>> Sequential version.\n" << std::endl;
#endif

  // Set up timers
  FTic time;
  FTic counter;

  // Print out some information on the input file and octree setup
  outputFileOctreeInfo(TreeHeight, SubTreeHeight, aboveTree, NbThreads, filename, periodicCondition);

  //! Currently calling loader twice, first time to populate octree, second time to add momenta. Reason: Right now you need populated octree to call container particle class which you need to add momenta.
  FFmaGenericLoader<FReal> loader(filename);
  geometryClass DevGeometry(&loader);
  
  OctreeClass tree(TreeHeight, SubTreeHeight, DevGeometry.BoxWidth, DevGeometry.BoxCenter);
  // Read particles and insert them in octree
  time.tic();
  //loadParticlesInOctree(&loader, TreeHeight, SubTreeHeight, &tree);
  time.tac();
  std::cout << "Done  " << "(@Creating and Inserting Particles = " << time.elapsed() << " s) ." << std::endl;


  EMCGenericLoader<FReal> testloader(filename);
  FSize Nb = testloader.getNumberOfParticles();
  FPoint<FReal> momenta;
  FReal *TotalMomenta;
  TotalMomenta = new FReal[3*Nb] ;
  LoadParticles(&testloader, TreeHeight, SubTreeHeight, &tree, &momenta, TotalMomenta);

  //* initialize meterial /scattering parameters

  mat_paramClass mat_par; // initialization of mat_par class,
  scat_paramClass scat_par(0,1); // Initialization of scattering after material parameters !! Acoustic, polar scattering 
  
  ScatteringTable(&scat_par, &mat_par);


  //* ------------- kspace and real space initalization --------------

  OctreeClass::Iterator octreeIterator(&tree);
  octreeIterator.gotoBottomLeft();

  int i, j, N;
  N = loader.getNumberOfParticles(); // total number of particles.
  i = 0; 
  {
  do{  
    ContainerClass* particles = octreeIterator.getCurrentListTargets();

    for(FSize idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart){ 


      k_inreader(particles, TotalMomenta,  &scat_par, &mat_par, idxPart, i);  // Change momenta

 //     if (particles->getPhysicalValues()[idxPart] > 0) // Charge index for particle indentifier
//        j++;

    //  Particle_Identifier(particles, idxPart, j, N); // give index to particles to determine their role (carrier or ion)


      i++;
      }
  
  } while(octreeIterator.moveRight());
  }  
  
  delete[] TotalMomenta;

  //* ----------------------------------------------------------------


  std::cout << "box dimension: " << DevGeometry.get_boxdim() << "\n \n" << std::endl;

  double pulsesize = 1.0;

  int totaltime = 10001;
  int tpulse = totaltime/10;
  int twrite = totaltime/200; //totaltime / 100;
  for(int tsim = 0 ; tsim < totaltime; tsim ++) // total simulation time in fs
  {
    std::cout << "Start of " << tsim << "'th FMM run" << std::endl;
    ////////////////////////////////////////////////////////////////////
    //
    //    Execute FMM Algorithm
    //
    ////////////////////////////////////////////////////////////////////

 // { // -----------------------------------------------------
    std::cout << "\n" << interpolationType << "  FMM (ORDER= "<< ORDER << ") ... " << std::endl;
    
    // reset all particle forces
    resetParticleLocation(&tree);
    
    octreeIterator.gotoBottomLeft();     

    const MatrixKernelClass  MatrixKernel;
    time.tic();

    std::unique_ptr<KernelClass> kernelsNoPer(new KernelClass(TreeHeight, DevGeometry.BoxWidth, DevGeometry.BoxCenter,&MatrixKernel));
    //
    FmmClass algoNoPer(&tree, kernelsNoPer.get());
    // periodic FMM algorithm
    FmmClassPer algoPer(&tree, aboveTree);
    KernelClass kernelsPer(algoPer.extendedTreeHeight(), algoPer.extendedBoxWidth(),
                           algoPer.extendedBoxCenter(),&MatrixKernel);
    algoPer.setKernel(&kernelsPer);
    //
    FAbstractAlgorithm * algorithm  = nullptr;
    FAlgorithmTimers   * timer      = nullptr;
    if(! periodicCondition) { // Non periodic case
        algorithm  = &algoNoPer ;
        timer      = &algoNoPer ;
        std::cout << "\n periodic is off \n" ;
        
      }
    else {                    // Periodic case
        algorithm  = &algoPer ;
        timer      = &algoPer ;
        std::cout << "\n periodic is on \n" ;
      }


    algorithm->execute();   // Here the call of the FMM algorithm

    time.tac();
    std::cout << "Timers Far Field \n"
              << "P2M " << timer->getTime(FAlgorithmTimers::P2MTimer) << " seconds\n"
              << "M2M " << timer->getTime(FAlgorithmTimers::M2MTimer) << " seconds\n"
              << "M2L " << timer->getTime(FAlgorithmTimers::M2LTimer) << " seconds\n"
              << "L2L " << timer->getTime(FAlgorithmTimers::L2LTimer) << " seconds\n"
              << "P2P and L2P " << timer->getTime(FAlgorithmTimers::NearTimer) << " seconds\n"
	      << std::endl;

    std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;




  //}
  // -----------------------------------------------------
  std::cout << " Particle positions updated according to new field\n" << std::endl;
  // Arrange Section--------------------------------------------------------------------------------*//
  
   



        std::cout << "Working on particles ..." << std::endl;
        counter.tic();

        {
          
         //   OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{  

                ContainerClass* particles2 = octreeIterator.getCurrentListTargets();
                particles2 = octreeIterator.getCurrentListTargets();
                // Advance each particle
               //----------------------------- MONTE-CARLO PROCEDURE ----------------------------------------------------------

                 

                for(FSize idxPart = 0; idxPart < particles2->getNbParticles() ; ++idxPart){ 

		//  if (particles2->getPhysicalValues(1,4)[idxPart] != 3) 
                  free_flight_scatter(DevGeometry, particles2, &scat_par, &mat_par, idxPart, tsim, tpulse); 
                  //std::cout << " \n \n particle potential before: " << *particles2->getPotentials();

                } // end for Advance each particle
            } while(octreeIterator.moveRight());
        }
        counter.tac();
        std::cout << "Done  " << "(@Moving = " << counter.elapsed() << "s)." << std::endl;

        const long repetitions = algoPer.theoricalRepetition();
        const long totalRepeatedBox = repetitions * repetitions * repetitions;
        std::cout << "The box is repeated " << repetitions << " there are " << totalRepeatedBox << " boxes in total\n";
        const long long NbParticlesEntireSystem = loader.getNumberOfParticles() * totalRepeatedBox;
        std::cout << "The total number of particles is "  << NbParticlesEntireSystem << "\n";
        


        // -- Pulse -- If particle is hole or electron it will be excited, until the pulsesize is reached.  Ions are excluced by if-statement. 
        
        if (tsim == tpulse ){
          octreeIterator.gotoBottomLeft();
          i = 0;
          {
          do{ 
          
          ContainerClass* particles3 = octreeIterator.getCurrentListTargets();

          for(FSize idxPart = 0; idxPart <  particles3->getNbParticles() ; ++idxPart)
          {
	   // if (particles3->getPhysicalValues(1,4)[idxPart] != 3)	
           // {
	      deltaE(particles3, DevGeometry, &scat_par, &mat_par, idxPart); 
	      i++;
          //  } 
          }
            
          }while(octreeIterator.moveRight() && i < pulsesize * N);
          }
        std::cout << "\n Pulse initiated. \n" << std::endl;
        }
        // ----- Writer ----------------
        
        
        
        if(tsim % twrite == 0){      
          
         std::string nombre;

         nombre = "../../../../data/p274598/Output/rr_noPI-10001-1e16_1e7-0.2-10-1.0/" + std::to_string(tsim) + ".fma";
         EMCGenericWriter<FReal> writeur(nombre);
         writeur.writeDataFromOctree(&tree, loader.getNumberOfParticles());
         //Input_Output.output_writer(&tree, tsim);
                
        }      

        
        // ---------------------Arranging -----------------------
        if(! periodicCondition) { // Non periodic case
          ArrangerClass arrange(&tree);
          std::cout << "Arrange ..." << std::endl;
          counter.tic();
          arrange.rearrange();
          counter.tac();
          std::cout << "Done  " << "(@Arrange = " << counter.elapsed() << "s)." << std::endl;

        
        }
        else {                    // Periodic case

          ArrangerClassPeriodic Arrange(&tree);
          std::cout << "Arrange periodic..." << std::endl;
          counter.tic();
          Arrange.rearrange();
          counter.tac();
          std::cout << "Done  " << "(@Arrange = " << counter.elapsed() << "s)." << std::endl;
        }
      
        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////



        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////


  }

  // ------------------------------------------
  // Write output using EMCGenericWriter class
  // ------------------------------------------

  std::cout <<std::endl<<"End of run " << std::endl<<std::endl;

  if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
    std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.fma"));
    EMCGenericWriter<FReal> writer(name) ;
    writer.writeDataFromOctree(&tree,loader.getNumberOfParticles());
  }

  //------------------------------------------------------------------------------------------------------

  
  return 0;
}

