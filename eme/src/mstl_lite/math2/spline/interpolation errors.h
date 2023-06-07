
#include "mstl.h"

// *******************************************************************
// *******************************************************************
// *******************************************************************

#if !defined(_SPLINE_ERRORS)
#define _SPLINE_ERRORS

namespace interpolation{

class NO_FIRST_SOLUTION;

class NO_TH_SOLUTION;

class NO_MORE_SOLUTIONS;

class NO_OVERLAP;

class TOO_FEW_POINTS;

// *******************************************************************
// *******************************************************************
// *******************************************************************

class NO_FIRST_SOLUTION : public std::runtime_error, public BASIC_ERROR
{ public : NO_FIRST_SOLUTION(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : runtime_error(std::string("No first solution")), BASIC_ERROR(WHICH,WHO,WHERE){;}
           NO_FIRST_SOLUTION(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : runtime_error(std::string("No first solution")), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

class NO_TH_SOLUTION : public std::runtime_error, public BASIC_ERROR
{ public : NO_TH_SOLUTION(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : runtime_error(std::string("No TH solution")), BASIC_ERROR(WHICH,WHO,WHERE){;}
           NO_TH_SOLUTION(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : runtime_error(std::string("No TH solution")), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

class NO_MORE_SOLUTIONS : public std::runtime_error, public BASIC_ERROR
{ public : NO_MORE_SOLUTIONS(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : runtime_error(std::string("No more solutions")), BASIC_ERROR(WHICH,WHO,WHERE){;}
           NO_MORE_SOLUTIONS(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : runtime_error(std::string("No more solutions")), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

class NO_OVERLAP : public std::runtime_error, public BASIC_ERROR
{ public : NO_OVERLAP(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : runtime_error(std::string("No overlap")), BASIC_ERROR(WHICH,WHO,WHERE){;}
           NO_OVERLAP(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : runtime_error(std::string("No overlap")), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

class TOO_FEW_POINTS : public std::invalid_argument, public BASIC_ERROR
{ public : TOO_FEW_POINTS(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : invalid_argument(std::string("Too few points")), BASIC_ERROR(WHICH,WHO,WHERE){;}
           TOO_FEW_POINTS(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : invalid_argument(std::string("Too few points")), BASIC_ERROR(WHICH,WHO,WHERE){;}
};

} // end of namespace

#endif
