#if !defined(_INTERPOLATION_DATA)
#include "interpolation data.h"
#endif

#if !defined(_INTERPOLATION_DATA_MEMBER_TEMPLATES)
#define _INTERPOLATION_DATA_MEMBER_TEMPLATES

namespace interpolation{

template<typename UNARYFUNCTOR> 
void XDATA_SINGLESET::TransformX(const UNARYFUNCTOR &UF)
     { 
       #if !defined(DNDEBUG)
       if(Empty()==true){ throw EMPTY("TransformX",CLASS_NAME,14);}       
       #endif
       std::vector<double> XT(NX()); 
       try{ int i;
            for(i=1;i<=static_cast<int>(NX());i++){ XT[i-1]=UF(_pX(i));} 
           }
       catch(...){ throw FUNCTION_ERROR("TransformX",CLASS_NAME,20);}
       _pSetX(XT);
      }

// *******************
// *******************
// *******************

template<typename UNARYFUNCTOR> 
void XDATA_MULTIPLESETS::TransformX(int i,const UNARYFUNCTOR &UF)
     { 
       #if !defined(DNDEBUG)
       if(Empty()==true){ throw EMPTY("TransformX",CLASS_NAME,32);}
       if(Empty(i)==true){ throw EMPTY("TransformX",CLASS_NAME,33);}
       #endif
       std::vector<double> XT(NX(i));
       try{ int j;
            for(j=1;j<=static_cast<int>(NX(i));j++){ XT[j-1]=UF(_pX(i,j));} 
           }
       catch(...){ throw FUNCTION_ERROR("TransformX",CLASS_NAME,39);}
       _pSetX(i,XT);
      }

template<typename UNARYFUNCTOR> 
void XDATA_MULTIPLESETS::TransformX(const UNARYFUNCTOR &UF)
     { 
       #if !defined(DNDEBUG)
       if(Empty()==true){ throw EMPTY("TransformX",CLASS_NAME,45);}
       #endif
       for(int i=1;i<=static_cast<int>(NXSets());i++)
              { 
                #if !defined(DNDEBUG)
                if(Empty(i)==true){ throw EMPTY("TransformX",CLASS_NAME,52);}
                #endif
                std::vector<double> XT(NX(i));
                try{ int j;
                     for(j=1;j<=static_cast<int>(NX(i));j++){ XT[j-1]=UF(_pX(i,j));} 
                    }
                catch(...){ throw FUNCTION_ERROR("TransformX",CLASS_NAME,7);}
                _pSetX(i,XT);
              }
      }

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<size_t ND> std::string YDATA_SINGLESET<ND>::CLASS_NAME("YDATA_SINGLESET<ND>");
//std::string YDATA_SINGLESET<1>::CLASS_NAME("YDATA_SINGLESET<1>");

template<size_t ND> 
double YDATA_SINGLESET<ND>::_pY(std::vector<int> i) const
                 { for(int a=0;a<=static_cast<int>(i.size())-1;a++){ i[a]--;}
                   return y[i];
                  }

template<size_t ND> 
void YDATA_SINGLESET<ND>::_pSetY(std::vector<int> i,double Y) const
                 { for(int a=0;a<=static_cast<int>(i.size())-1;a++){ i[a]--;}
                   y[i]=Y;
                  }

// *******************************************************

template<size_t ND>
template<typename UNARYFUNCTOR> 
void YDATA_SINGLESET<ND>::_pTransformY(std::vector<int> i,int j,const UNARYFUNCTOR &UF) 
     { if(j==static_cast<int>(ND)+1){ _pSetY(i,UF(_pY(i)));} 
       else{ for(int a=1;a<=static_cast<int>(ND);a++){ i[j-1]=a; _pTransformY(i,j+1,UF);} }
      }

// shift, rescale the data      
template<size_t ND> 
void YDATA_SINGLESET<ND>::_pShiftY(std::vector<int> i,int j,const double &S)
     { if(j==static_cast<int>(ND)+1){ _pSetY(i,_pY(i)+S);}
       else{ for(int a=1;a<=static_cast<int>(NY(j));a++){ i[j-1]=a; _pShiftY(i,j+1,S);} }
      }

template<size_t ND> 
void YDATA_SINGLESET<ND>::_pRescaleY(std::vector<int> i,int j,const double &S)
     { if(j==static_cast<int>(ND)+1){ _pSetY(i,_pY(i)*S);}
       else{ for(int a=1;a<=static_cast<int>(NY(j));a++){ i[j-1]=a; _pRescaleY(i,j+1,S);} }
      }     

template<size_t ND> 
void YDATA_SINGLESET<ND>::_pAbsY(std::vector<int> i,int j)
     { if(j==static_cast<int>(ND)+1){ _pSetY(i,Sign(_pY(i))*_pY(i));}
       else{ for(int a=1;a<=static_cast<int>(NY(j));a++){ i[j-1]=a; _pAbsY(i,j+1);} }
      }

// *******************************************************

template<size_t ND> 
double YDATA_SINGLESET<ND>::Y(std::vector<int> i) const
       { 
         #if !defined(DNDEBUG)
         if(i.size()!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i.size()),ND,std::string("Y"),CLASS_NAME);}
         #endif       
         std::vector<size_t> D=NY(); 
         #if !defined(DNDEBUG)         
         for(int a=0;a<=static_cast<int>(ND)-1;a++){ if(i[a]<1 || i[a]>D[a]){ throw OUT_OF_RANGE<int>(i[a],1,D[a],std::string("Y"),CLASS_NAME);} }
         #endif
         return _pY(i);         
        }

template<size_t ND> 
void YDATA_SINGLESET<ND>::SetY(std::vector<int> i,double Y) const
       { 
         #if !defined(DNDEBUG)
         if(i.size()!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i.size()),ND,std::string("Y"),CLASS_NAME);}
         #endif
         std::vector<size_t> D=NY(); 
         #if !defined(DNDEBUG)         
         for(int a=0;a<=static_cast<int>(ND)-1;a++){ if(i[a]<1 || i[a]>D[a]){ throw OUT_OF_RANGE<int>(i[a],1,D[a],std::string("SetY"),CLASS_NAME);} }
         #endif
         _pSetY(i,Y);
        }

// *******************************************************

template<size_t ND> template<typename UNARYFUNCTOR> 
void YDATA_SINGLESET<ND>::TransformY(const UNARYFUNCTOR &UF)
     { 
       #if !defined(DNDEBUG)
       if(Empty()==true){ throw EMPTY("TransformY",CLASS_NAME,9);}       
       #endif
       std::vector<int> i(ND);
       int j=1;
       for(int a=1;a<=static_cast<int>(NY(j));a++){ i[j-1]=a; _pTransformY(i,j+1,UF);}
      }

template<size_t ND> 
void YDATA_SINGLESET<ND>::ShiftY(const double &S)
        { 
          #if !defined(DNDEBUG)
          if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME);}
          #endif
          std::vector<int> i(ND);
          int j=1;
          for(int a=1;a<=static_cast<int>(NY(j));a++){ i[j-1]=a; _pShiftY(i,j+1,S);} 
         }

template<size_t ND> 
void YDATA_SINGLESET<ND>::RescaleY(const double &S)
        { 
          #if !defined(DNDEBUG)
          if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME);}
          #endif
          std::vector<int> i(ND);
          int j=1;
          for(int a=1;a<=static_cast<int>(NY(j));a++){ i[j-1]=a; _pRescaleY(i,j+1,S);} 
         }

template<size_t ND> 
void YDATA_SINGLESET<ND>::AbsY(void)
        { 
          #if !defined(DNDEBUG)
          if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME);}
          #endif
          std::vector<int> i(ND);
          int j=1;
          for(int a=1;a<=static_cast<int>(NY(j));a++){ i[j-1]=a; _pAbsY(i,j+1);} 
         }

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<typename UNARYFUNCTOR> 
void YDATA_SINGLESET<1>::TransformY(const UNARYFUNCTOR &UF)
     { 
       #if !defined(DNDEBUG)
       if(Empty()==true){ throw EMPTY("TransformY",CLASS_NAME,9);}       
       #endif
       std::vector<double> YT(NY()); 
       try{ int i;
            for(i=1;i<=static_cast<int>(NY());i++){ YT[i-1]=UF(_pY(i));} 
           }
       catch(...){ throw FUNCTION_ERROR("TransformY",CLASS_NAME,7);}
       _pSetY(YT);
      }

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<typename UNARYFUNCTOR> 
void YDATA_MULTIPLESETS::TransformY(int i,const UNARYFUNCTOR &UF) 
     { 
       #if !defined(DNDEBUG)
       if(Empty()==true){ throw EMPTY("TransformY",CLASS_NAME,9);}
       if(Empty(i)==true){ throw EMPTY("TransformY",CLASS_NAME,9);}
       #endif
       std::vector<double> YT(NY(i));
       try{ int j;
            for(j=1;j<=static_cast<int>(NY(i));j++){ YT[j-1]=UF(_pY(i,j));} 
           }
       catch(...){ throw FUNCTION_ERROR("TransformY",CLASS_NAME,7);}
       _pSetY(i,YT);
      }

template<typename UNARYFUNCTOR> 
void YDATA_MULTIPLESETS::TransformY(const UNARYFUNCTOR &UF) 
     { 
       #if !defined(DNDEBUG)
       if(Empty()==true){ throw EMPTY("TransformY",CLASS_NAME,9);}
       #endif
       for(int i=1;i<=static_cast<int>(NYSets());i++)
          { 
            #if !defined(DNDEBUG)
            if(Empty(i)==true){ throw EMPTY("TransformY",CLASS_NAME,9);}
            #endif
            std::vector<double> YT(NY(i));
            try{ int j;
                 for(j=1;j<=static_cast<int>(NY(i));j++){ YT[j-1]=UF(_pY(i,j));} 
                }
            catch(...){ throw FUNCTION_ERROR("TransformY",CLASS_NAME,7);}
            _pSetY(i,YT);
          }
      }

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<size_t ND> std::string SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::CLASS_NAME("SPLINE_MULTIPLEXSETS_SINGLEYSET");

template<size_t ND> 
double SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::_pA(std::vector<int> i,std::vector<int> j) const
       { for(int x=0;x<=static_cast<int>(ND)-1;x++){ i[x]--;}
         return a[i][j];
        }

template<size_t ND> 
TARRAY<double,ND> SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::_pA(std::vector<int> i) const
       { for(int x=0;x<=static_cast<int>(ND)-1;x++){ i[x]--;}
         return a[i];
        }

template<size_t ND> 
void SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::_pSetA(std::vector<int> i,std::vector<int> j,double D) const
       { for(int x=0;x<=static_cast<int>(ND)-1;x++){ i[x]--;}
         a[i][j]=D;
         SetFitted(false);
        }

template<size_t ND> 
void SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::_pSetA(std::vector<int> i,TARRAY<double,ND> D) const
       { for(int x=0;x<=static_cast<int>(ND)-1;x++){ i[x]--;}
         a[i]=D;
         SetFitted(false);
        }

// *******************************

template<size_t ND> 
void SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::CreateASet(std::vector<int> i,std::vector<size_t> NP)
     { for(int x=0;x<=static_cast<int>(ND)-1;x++){ i[x]--;}
       a[i]=TARRAY<double,ND>(NP);
       SetFitted(false);
      }

template<size_t ND> 
void SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::DestroyA(void) const
     { std::vector<size_t> N(ND,0);
       a=TARRAY<TARRAY<double,ND>,ND>(N,TARRAY<double,ND>(N));
       SetFitted(false);
      }

// *******************************

template<size_t ND> 
std::vector<size_t> SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::NParameters(std::vector<int> i) const
     { for(int x=0;x<=static_cast<int>(i.size())-1;x++){ i[x]--;}
       return a[i].Dimenions();
      }

template<size_t ND> 
double SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::A(std::vector<int> i,std::vector<int> j) const
       { 
         #if !defined(DNDEBUG)
         if(i.size()!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i.size()),ND,std::string("A"),CLASS_NAME);}
         if(static_cast<int>(j.size())!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(j.size()),ND,std::string("A"),CLASS_NAME);}
         #endif
         std::vector<size_t> D=NA(); 
         #if !defined(DNDEBUG)          
         for(int a=0;a<=static_cast<int>(ND)-1;a++){ if(i[a]<1 || i[a]>D[a]){ throw OUT_OF_RANGE<int>(i[a],1,D[a],std::string("A"),CLASS_NAME);} }
         #endif
         return _pA(i,j);
        }

template<size_t ND> 
TARRAY<double,ND> SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::A(std::vector<int> i) const
       { 
         #if !defined(DNDEBUG)
         if(i.size()!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i.size()),ND,std::string("A"),CLASS_NAME);}
         #endif
         std::vector<size_t> D=NA(); 
         #if !defined(DNDEBUG)
         for(int a=0;a<=static_cast<int>(ND)-1;a++){ if(i[a]<1 || i[a]>D[a]){ throw OUT_OF_RANGE<int>(i[a],1,D[a],std::string("A"),CLASS_NAME);} }
         #endif 
         return _pA(i);
        }

template<size_t ND> 
void SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::SetA(std::vector<int> i,std::vector<int> j,double A) const
       { 
         #if !defined(DNDEBUG)
         if(i.size()!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i.size()),ND,std::string("SetA"),CLASS_NAME);}
         if(static_cast<int>(j.size())!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(j.size()),ND,std::string("SetA"),CLASS_NAME);}
         #endif
         std::vector<size_t> D=NA(); 
         for(int a=0;a<=static_cast<int>(ND)-1;a++){ if(i[a]<1 || i[a]>D[a]){ throw OUT_OF_RANGE<int>(i[a],1,D[a],std::string("SetA"),CLASS_NAME);} }
         _pSetA(i,j,A);
        }

template<size_t ND> 
void SPLINE_MULTIPLEXSETS_SINGLEYSET<ND>::SetA(std::vector<int> i,TARRAY<double,ND> A) const
       { 
         #if !defined(DNDEBUG)
         if(i.size()!=ND){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i.size()),ND,std::string("SetA"),CLASS_NAME);}
         #endif
         std::vector<size_t> D=NA(); 
         #if !defined(DNDEBUG)        
         for(int a=0;a<=static_cast<int>(ND)-1;a++){ if(i[a]<1 || i[a]>D[a]){ throw OUT_OF_RANGE<int>(i[a],1,D[a],std::string("SetA"),CLASS_NAME);} }
         #endif
         _pSetA(i,A);
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
                #if !defined(DNDEBUG)
                if(D==0.){ throw DIVISION_BY_ZERO("/=(UNARYMEMBERFUNCTION)");}
                #endif 
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
                #if !defined(DNDEBUG)
                if(D==0.){ throw DIVISION_BY_ZERO("/=(BINARYMEMBERFUNCTION)");}
                #endif
                _pSetY(i,_pY(i)/D);
               }
           return (*this);
          }        
*/
} // end of namespace interpolation

#endif

