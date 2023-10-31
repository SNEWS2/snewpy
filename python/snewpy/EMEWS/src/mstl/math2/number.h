
#include <cmath>

#include "mstl.h"

#if !defined(_NUMBERS)
#define _NUMBERS

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

// These constants are defined in the math.h header
// M_E, M_LOG2E, M_LOG10E, M_LN2, M_LN10, 
// M_PI, M_PI_2, M_PI_4, M_1_PI, M_2_PI, M_2_SQRTPI
// M_SQRT2, M_SQRT1_2

#define M_SQRTPI 1.772453850905516027298167      // sqrt of pi
#define M_1_SQRTPI 0.564189583547756286948       // reciprocal of sqrt pi 
#define M_2PI 6.283185307179586476925286766559   // 2 pi
#define M_4PI 2.*M_2PI                           // 4 pi
#define M_PI2 M_PI*M_PI                          // pi ^2
#define M_PI3 M_PI*M_PI2                         // pi ^3
#define M_PI4 M_PI2*M_PI2                        // pi ^4 

#define M_PI_3 M_PI / 3.        // pi / 3
#define M_2PI_3 M_2PI / 3.      // 2pi / 3
#define M_4PI_3 2.* M_2PI_3     // 4pi / 3

#define M_SQRT3 1.7320508075688772935  // sqrt of 3
#define M_SQRT5 2.2360679774997896964  // sqrt of 5
#define M_SQRT6 M_SQRT2*M_SQRT3        // sqrt of 6

#define M_CBRT2 1.259921049894873      // cube root of 2

#define M_COS0 1.0                           // cosine of 0 degrees
#define M_COS15 (M_SQRT6 + M_SQRT2) / 4.     // cosine of 15 degrees
#define M_COS30 M_SQRT3 / 2.                 // cosine of 30 degrees
#define M_COS45 M_SQRT1_2                    // cosine of 45 degrees
#define M_COS60 0.5                          // cosine of 60 degrees
#define M_COS75 (M_SQRT6 - M_SQRT2) / 4.     // cosine of 75 degrees
#define M_COS90 0.0                          // cosine of 90 degrees
#define M_COS105 -(M_SQRT6 - M_SQRT2) / 4.   // cosine of 105 degrees
#define M_COS120 -0.5                        // cosine of 120 degrees
#define M_COS135 -M_SQRT1_2                  // cosine of 135 degrees
#define M_COS150 -M_SQRT3 / 2.               // cosine of 150 degrees
#define M_COS165 -(M_SQRT6 + M_SQRT2) / 4.   // cosine of 165 degrees
#define M_COS180 -1.0                        // cosine of 180 degrees
#define M_COS195 -(M_SQRT6 + M_SQRT2) / 4.   // cosine of 195 degrees
#define M_COS210 -M_SQRT3 / 2.               // cosine of 210 degrees
#define M_COS225 -M_SQRT1_2                  // cosine of 225 degrees
#define M_COS240 -0.5                        // cosine of 240 degrees
#define M_COS255 -(M_SQRT6 - M_SQRT2) / 4.   // cosine of 255 degrees
#define M_COS270 0.0                         // cosine of 270 degrees
#define M_COS285 (M_SQRT6 - M_SQRT2) / 4.    // cosine of 285 degrees
#define M_COS300 0.5                         // cosine of 300 degrees
#define M_COS315 M_SQRT1_2                   // cosine of 215 degrees
#define M_COS330 M_SQRT3 / 2.                // cosine of 330 degrees
#define M_COS345 (M_SQRT6 + M_SQRT2) / 4.    // cosine of 345 degrees
#define M_COS360 1.0                         // cosine of 360 degrees

#define M_SIN0 0.0                           // sine of 0 degrees
#define M_SIN15 (M_SQRT6 - M_SQRT2) / 4.     // sine of 15 degrees
#define M_SIN30 0.5                          // sine of 30 degrees
#define M_SIN45 M_SQRT1_2                    // sine of 45 degrees
#define M_SIN60 M_SQRT3 / 2.                 // sine of 60 degrees
#define M_SIN75 (M_SQRT6 + M_SQRT2) / 4.     // sine of 75 degrees
#define M_SIN90 1.0                          // sine of 90 degrees
#define M_SIN105 (M_SQRT6 + M_SQRT2) / 4.    // sine of 105 degrees
#define M_SIN120 M_SQRT3 / 2.                // sine of 120 degrees
#define M_SIN135 M_SQRT1_2                   // sine of 135 degrees
#define M_SIN150 0.5                         // sine of 150 degrees
#define M_SIN165 (M_SQRT6 - M_SQRT2) / 4.    // sine of 165 degrees
#define M_SIN180 0.0                         // sine of 180 degrees
#define M_SIN195 -(M_SQRT6 - M_SQRT2) / 4.   // sine of 195 degrees
#define M_SIN210 -0.5                        // sine of 210 degrees
#define M_SIN225 -M_SQRT1_2                  // sine of 225 degrees
#define M_SIN240 -M_SQRT3 / 2.               // sine of 240 degrees
#define M_SIN255 -(M_SQRT6 + M_SQRT2) / 4.   // sine of 255 degrees
#define M_SIN270 -1.0                        // sine of 270 degrees
#define M_SIN285 -(M_SQRT6 + M_SQRT2) / 4.   // sine of 285 degrees
#define M_SIN300 -M_SQRT3 / 2.               // sine of 300 degrees
#define M_SIN315 -M_SQRT1_2                  // sine of 215 degrees
#define M_SIN330 -0.5                        // sine of 330 degrees
#define M_SIN345 -(M_SQRT6 - M_SQRT2) / 4.   // sine of 345 degrees
#define M_SIN360 0.0                         // sine of 360 degrees

#define M_GAMMA1 0
#define M_GAMMA2 1
#define M_GAMMA3 2
#define M_GAMMA4 6
#define M_GAMMA5 24
#define M_GAMMA6 120
#define M_GAMMA7 720
#define M_GAMMA8 5040
#define M_GAMMA9 40320
#define M_GAMMA10 362880
#define M_GAMMA11 3628800 
#define M_GAMMA12 39916800
#define M_GAMMA13 479001600
#define M_GAMMA14 6227020800
#define M_GAMMA15 87178291200
#define M_GAMMA16 1307674368000
#define M_GAMMA17 20922789888000
#define M_GAMMA18 355687428096000
#define M_GAMMA19 6402373705728000
#define M_GAMMA20 121645100408832000

// pure numbers
const double EulerMascheroni=0.577215664901532861;
const double zeta2=M_PI2/6., zeta3=1.2020569031595942853997381, zeta4=M_PI4/90.;
const double Catalan=0.915965594177219015054603514932384110774;

#endif
