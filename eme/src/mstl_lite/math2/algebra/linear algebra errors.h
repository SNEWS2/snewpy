
#include "mstl.h"

// *******************************************************************
// *******************************************************************
// *******************************************************************

#if !defined(_LINALG_ERRORS)
#define _LINALG_ERRORS

// ********************* //
// LINEAR ALGEBRA ERRORS //
// ********************* //

class INCORRECT_FORM;

class NOT_SQUARE;

class NOT_SYMMETRIC;

class NOT_INVERTABLE;

class SINGULAR;

// *******************************************************************
// *******************************************************************
// *******************************************************************

class INCORRECT_FORM : public std::invalid_argument, public BASIC_ERROR
{ public : INCORRECT_FORM(const std::string &WHAT=std::string("Incorrect form"),const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(WHAT), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

class NOT_SQUARE : public std::invalid_argument, public BASIC_ERROR
{ public : NOT_SQUARE(const std::string &WHAT=std::string("Not square"),const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(WHAT), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

class NOT_SYMMETRIC : public std::invalid_argument, public BASIC_ERROR
{ public : NOT_SYMMETRIC(const std::string &WHAT=std::string("Not symmetric"),const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(WHAT), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

class NOT_INVERTABLE : public std::runtime_error, public BASIC_ERROR
{ public : NOT_INVERTABLE(const std::string &WHAT=std::string("Not invertable"),const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::runtime_error(WHAT), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

class SINGULAR : public std::runtime_error, public BASIC_ERROR
{ public : SINGULAR(const std::string &WHAT=std::string("Singular matrix"),const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::runtime_error(WHAT), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

#endif


