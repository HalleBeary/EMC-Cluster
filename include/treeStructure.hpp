#ifndef treeStructure_HPP
#define treeStructure_HPP

#include <string> 
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"





#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel_Covariance.hpp"

//
template<typename FReal, int ORDER> 
using FInterpolationCell = FChebCell<FReal, ORDER>; // To change interpolation type, change FInterpolationCell and FInterpolationKernel
template<typename FReal, typename GroupCellClass,
	 typename GroupContainerClass,
	 typename MatrixKernelClass, int ORDER>  
using FInterpolationKernel =  FChebSymKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> ;
							 //FChebSymKernel<FReal,GroupCellClass, GroupContainerClass, MatrixKernelClass, ORDER> ;
const std::string interpolationType("Uniform interpolation");

#include <iostream>
#include <iomanip>
#include <memory>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"
#include "Utils/FTic.hpp" // copy
#include "Utils/FPoint.hpp" // copy

#include "Files/FFmaGenericLoader.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Components/FSimpleLeaf.hpp"

/// definition of the common tree structure
static constexpr unsigned ORDER = 2;
static constexpr unsigned core_width = 0.1;

using FReal = double;
using CellClass      = FInterpolationCell<FReal, ORDER>;
using CellUpClass    = typename CellClass::multipole_t;
using ContainerClass = FP2PParticleContainerIndexed<FReal,6>;
using LeafClass      = FSimpleLeaf<FReal,  ContainerClass >   ;
using OctreeClass    = FOctree<FReal, CellClass,ContainerClass,LeafClass>;

#include "Components/FBasicParticleContainer.hpp"
#include "Components/FBasicCell.hpp"

// Arranger
#include "Arranger/FBasicParticleContainerIndexedMover.hpp"
#include "Arranger/FOctreeArranger.hpp"
#include "Arranger/FArrangerPeriodic.hpp"

// Octree
#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"
//
#include "Kernels/Interpolation/FCutOffKernel.hpp"

using MatrixKernelClass = FInterpMatrixKernelAPLUSR<FReal> ;   // Set the kernel you want
using KernelClass       = FInterpolationKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> ;

#ifdef _OPENMP
  #include "Core/FFmmAlgorithmPeriodic.hpp"
  #include "Core/FFmmAlgorithmSectionTask.hpp"
#else
  #include "Core/FFmmAlgorithm.hpp"
#endif


#ifdef _OPENMP
using FmmClass    = FFmmAlgorithmSectionTask<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> ;
using FmmClassPer = FFmmAlgorithmPeriodic<FReal,OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;
#else
using FmmClass  = FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;
#endif

#endif //treeStructure_HPP
