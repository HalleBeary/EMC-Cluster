///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INASTEMPCOMPILECONFIG_H
#define INASTEMPCOMPILECONFIG_H

#ifdef INASTEMP_USE_STATIC_CONFIG
#error "It is forbiden to include InastempCompileConfig in case of static config, please protect inclusion with INASTEMP_USE_STATIC_CONFIG."
#endif

// Define all macros (ADD-NEW-HERE)
#define INASTEMP_USE_SCALAR

#define INASTEMP_USE_SSE3

#define INASTEMP_USE_SSSE3

#define INASTEMP_USE_SSE41

#define INASTEMP_USE_SSE42

#define INASTEMP_USE_AVX

#define INASTEMP_USE_AVX2

/* #undef INASTEMP_USE_AVX512COMMON */

/* #undef INASTEMP_USE_AVX512KNL */

/* #undef INASTEMP_USE_AVX512SKL */

/* #undef INASTEMP_USE_ALTIVEC */
/* #undef INASTEMP_USE_XL */

// Inform about best one
#define INASTEMP_AVX2_IS_BEST
#define INASTEMP_BEST_TYPE AVX2

#ifndef INASTEMP_NO_BEST_INCLUDE
#include "AVX2/InaVecAVX2Float.hpp"
#include "AVX2/InaVecAVX2Double.hpp"
#else
template <class RealType>
class InaVecAVX2;
#endif

template <class RealType>
using InaVecBestType = InaVecAVX2<RealType>;

using InaVecBestTypeFloat = InaVecAVX2<float>;
using InaVecBestTypeDouble = InaVecAVX2<double>;


#endif
