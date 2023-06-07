#if !defined(_DISCONTINUOUS)
#include "discontinuous.h"
#endif

#if !defined(_DISCONTINUOUS_MEMBERS)
#define _DISCONTINUOUS_MEMBERS

namespace interpolation{

template <typename IType> DISCONTINUOUS::DISCONTINUOUS(EXPRESSION<IType> const &C)
         : XDATA_SINGLESET(C.N()), YDATA_SINGLESET<1>(C.N())
         { SetFoundDomains(false);
           int i;
           for(i=1;i<=C.N();i++){ _pSetX(i,C.X(i)); _pSetY(i,C(X(i)));}
           SetFoundDomains(false);
          }


template <typename IType> DISCONTINUOUS& DISCONTINUOUS::operator=(EXPRESSION<IType> const &C)
         { CreateX(C.N()); CreateY(C.N());
           int i;
           for(i=1;i<=C.N();i++){ _pSetX(i,C.X(i)); _pSetY(i,C(X(i)));} 
           SetXDifferenced(false); 
           SetGraded(false); 
           SetFitted(false); 
           SetSorted(false); 
           SetYLimited(false); 
           SetXYLimited(false); 
           SetFoundDomains(false);
           return *this;
          }

} //end of namespace interpolation

#endif
