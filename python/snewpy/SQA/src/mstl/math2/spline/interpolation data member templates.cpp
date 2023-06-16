#if !defined(_INTERPOLATION_DATA_MEMBERS)
#include "interpolation data.h"
#endif

#if !defined(_DATA_INTERPOLATION_MEMBERS)
#define _DATA_INTERPOLATION_MEMBERS

namespace interpolation{

template<typename UNARYFUNCTOR> 
void XDATA_SINGLESET::TransformX(const UNARYFUNCTOR &UF)
     { if(Empty()==true){ throw EMPTY("TransformX",CLASS_NAME,12);}       
       std::vector<double> XT(NX()); 
       try{ int i;
            for(i=1;i<=NX();i++){ XT[i-1]=UF(_pX(i));} 
           }
       catch(...){ throw FUNCTION_ERROR("TransformX",CLASS_NAME,17);}
       _pSetX(XT);
      }

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<int ND> std::string YDATA_SINGLESET<ND>::CLASS_NAME("YDATA_SINGLESET");

template<int ND> 
double YDATA_SINGLESET<ND>::_pY(std::vector<int> i) const
                 { for(int a=0;a<=static_cast<int>(i.size())-1;a++){ i[a]--;}
                   return y[i];
                  }

template<int ND> 
void YDATA_SINGLESET<ND>::_pSetY(std::vector<int> i,double Y) const
                 { for(int a=0;a<=static_cast<int>(i.size())-1;a++){ i[a]--;}
                   y[i]=Y;
                  }

// *******************************************************

template<int ND>
template<typename UNARYFUNCTOR> 
void YDATA_SINGLESET<ND>::_pTransformY(std::vector<int> i,int j,const UNARYFUNCTOR &UF) 
     { if(j==ND+1){ _pSetY(i,UF(_pY(i)));} 
       else{ for(int a=1;a<=ND;a++){ i[j-1]=a; _pTransformY(i,j+1,UF);} }
      }

// shift, rescale the data      
template<int ND> 
void YDATA_SINGLESET<ND>::_pShiftY(std::vector<int> i,int j,const double &S)
     { if(j==ND+1){ _pSetY(i,_pY(i)+S);}
       else{ for(int a=1;a<=NY(j);a++){ i[j-1]=a; _pShiftY(i,j+1,S);} }
      }

template<int ND> 
void YDATA_SINGLESET<ND>::_pRescaleY(std::vector<int> i,int j,const double &S)
     { if(j==ND+1){ _pSetY(i,_pY(i)*S);}
       else{ for(int a=1;a<=NY(j);a++){ i[j-1]=a; _pRescaleY(i,j+1,S);} }
      }     

template<int ND> 
void YDATA_SINGLESET<ND>::_pAbsY(std::vector<int> i,int j)
     { if(j==ND+1){ _pSetY(i,Sign(_pY(i))*_pY(i));}
       else{ for(int a=1;a<=NY(j);a++){ i[j-1]=a; _pAbsY(i,j+1);} }
      }

// *******************************************************

template<int ND> 
double YDATA_SINGLESET<ND>::Y(std::vector<int> i) const
       { if(static_cast<int>(i.size())!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i.size()),ND,std::string("Y"),CLASS_NAME);}
         std::vector<int> D=NY(); 
         for(int a=0;a<=ND-1;a++){ if(i[a]<1 || i[a]>D[a]){ throw OUT_OF_RANGE<int>(i[a],1,D[a],std::string("Y"),CLASS_NAME);} }
         return _pY(i);         
        }

template<int ND> 
void YDATA_SINGLESET<ND>::SetY(std::vector<int> i,double Y) const
       { if(static_cast<int>(i.size())!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i.size()),ND,std::string("Y"),CLASS_NAME);}
         std::vector<int> D=NY(); 
         for(int a=0;a<=ND-1;a++){ if(i[a]<1 || i[a]>D[a]){ throw OUT_OF_RANGE<int>(i[a],1,D[a],std::string("SetY"),CLASS_NAME);} }
         _pSetY(i,Y);
        }

// *******************************************************

template<int ND> template<typename UNARYFUNCTOR> 
void YDATA_SINGLESET<ND>::TransformY(const UNARYFUNCTOR &UF)
     { if(Empty()==true){ throw EMPTY("TransformY",CLASS_NAME,9);}       
       std::vector<int> i(ND);
       int j=1;
       for(int a=1;a<=NY(j);a++){ i[j-1]=a; _pTransformY(i,j+1,UF);}
      }

template<int ND> 
void YDATA_SINGLESET<ND>::ShiftY(const double &S)
        { if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME);}
          std::vector<int> i(ND);
          int j=1;
          for(int a=1;a<=NY(j);a++){ i[j-1]=a; _pShiftY(i,j+1,S);} 
         }

template<int ND> 
void YDATA_SINGLESET<ND>::RescaleY(const double &S)
        { if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME);}
          std::vector<int> i(ND);
          int j=1;
          for(int a=1;a<=NY(j);a++){ i[j-1]=a; _pRescaleY(i,j+1,S);} 
         }

template<int ND> 
void YDATA_SINGLESET<ND>::AbsY(void)
        { if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME);}
          std::vector<int> i(ND);
          int j=1;
          for(int a=1;a<=NY(j);a++){ i[j-1]=a; _pAbsY(i,j+1);} 
         }

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<typename UNARYFUNCTOR> 
void YDATA_SINGLESET<1>::TransformY(const UNARYFUNCTOR &UF)
     { if(Empty()==true){ throw EMPTY("TransformY",CLASS_NAME,125);}       
       std::vector<double> YT(NY()); 
       try{ int i;
            for(i=1;i<=NY();i++){ YT[i-1]=UF(_pY(i));} 
           }
       catch(...){ throw FUNCTION_ERROR("TransformY",CLASS_NAME,130);}
       _pSetY(YT);
      }

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************
/*
template <typename UMFType> OPERATOR& OPERATOR::operator+=(const UNARYMEMBERFUNCTION<UMFType,double> &UMF)
        { int i;
                   for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)+UMF(_pX(i)) );} 
           return (*this);
          }

template <typename UMFType> OPERATOR& OPERATOR::operator-=(const UNARYMEMBERFUNCTION<UMFType,double> &UMF)
        { int i;
                   for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)-UMF(_pX(i)) );} 
           return (*this);
          }

template <typename UMFType> OPERATOR& OPERATOR::operator*=(const UNARYMEMBERFUNCTION<UMFType,double> &UMF)
        { int i;
                   for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)*UMF(_pX(i)) );} 
           return (*this);
          }

template <typename UMFType> OPERATOR& OPERATOR::operator/=(const UNARYMEMBERFUNCTION<UMFType,double> &UMF)
         { double D;
            int i;
            for(i=1;i<=NY();i++)
                 { D=UMF(_pX(i));
                    if(D==0.){ throw DIVISION_BY_ZERO("/=(UNARYMEMBERFUNCTION)");}
                    _pSetY(i,_pY(i)/D);
                  }
           return (*this);
          }  

template <typename BMFType> OPERATOR& OPERATOR::operator+=(const BINARYMEMBERFUNCTION<BMFType,double> &BMF)
        { int i;
                   for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)+BMF(_pX(i),_pY(i)) );} 
           return (*this);
         }

template <typename BMFType> OPERATOR& OPERATOR::operator-=(const BINARYMEMBERFUNCTION<BMFType,double> &BMF)
        { int i;
                   for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)-BMF(_pX(i),_pY(i)) );} 
           return (*this);
          }

template <typename BMFType> OPERATOR& OPERATOR::operator*=(const BINARYMEMBERFUNCTION<BMFType,double> &BMF)
        { int i;
                   for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)*BMF(_pX(i),_pY(i)) );} 
           return (*this);
          }

template <typename BMFType> OPERATOR& OPERATOR::operator/=(const BINARYMEMBERFUNCTION<BMFType,double> &BMF)
         { double D;
            int i;
            for(i=1;i<=NY();i++)
                   { D=BMF(_pX(i),_pY(i));
                      if(D==0.){ throw DIVISION_BY_ZERO("/=(BINARYMEMBERFUNCTION)");}
                      _pSetY(i,_pY(i)/D);
                    }
           return (*this);
          }        
*/
} // end of namespace interpolation

#endif

