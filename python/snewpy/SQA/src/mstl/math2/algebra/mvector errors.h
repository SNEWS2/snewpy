
#include "mstl.h"

// *******************************************************************
// *******************************************************************
// *******************************************************************

#if !defined(_MVECTOR_ERRORS)
#define _MVECTOR_ERRORS

// ***************//
// MVECTOR ERRORS //
// ***************//

class NOT_3D;

// ***************************************************

class NOT_3D : public BASIC_ERROR { public : NOT_3D(char *FUNCTION) : BASIC_ERROR("Not 3D vectors",FUNCTION,"MVECTOR") {;} };

#endif


