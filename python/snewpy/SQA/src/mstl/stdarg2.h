
#include<cstdarg>

#include "mstl.h"

#if !defined(_STDARG99)
#define _STDARG99

#ifdef _ICC
va_list va_make(std::va_list, ...);
#endif

#endif

