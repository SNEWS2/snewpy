
#include "stdarg2.h"

using std::va_list;

#ifdef _ICC
va_list va_make(va_list ap, ...){ va_start(ap,ap); return ap;}
#endif
