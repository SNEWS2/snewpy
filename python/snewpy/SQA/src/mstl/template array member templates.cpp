#if !defined(_TARRAY_MEMBERS)
#define _TARRAY_MEMBERS

// *******************************************************************
// *******************************************************************
// *******************************************************************

// rank >=3
template <class Type,size_t rank> template <size_t rank2>
void TARRAY<Type,rank>::Create(TARRAY<Type,rank2>*,std::vector<size_t> DIMENSIONS)
         { //assert(DIMENSIONS.size()==rank);
           #if !defined(DNDEBUG)
           if(DIMENSIONS.size()!=rank){ throw DIFFERENT_DIMENSIONS("Create(TARRAY<Type,rank2>*,std::vector<int>)","TARRAY<Type,rank>");}
           #endif
           dimensions=DIMENSIONS;

           if(dimensions[0]!=0)
             { empty=false;
               mt=new TARRAY<Type,rank-1>[dimensions[0]];
               DIMENSIONS.erase(DIMENSIONS.begin());

               for(int i=0;i<=static_cast<int>(dimensions[0])-1;i++){ mt[i].Create(&mt[i],DIMENSIONS);}
              }
           else{ empty=true; mt=NULL;}
          }

// partial specialization for rank 2 tensors
template <class Type,size_t rank> 
void TARRAY<Type,rank>::Create(TARRAY<Type,2>*,std::vector<size_t> DIMENSIONS)
         { //assert(DIMENSIONS.size()==2);
           #if !defined(DNDEBUG)
           if(DIMENSIONS.size()!=2){ throw DIFFERENT_DIMENSIONS("Create(TARRAY<Type,2>*,std::vector<int>)","TARRAY<Type,1>");}
           #endif
           dimensions=DIMENSIONS;

           if(dimensions[0]!=0)
             { empty=false;
               mt=new TARRAY<Type,1>[dimensions[0]];

               for(int i=0;i<=static_cast<int>(dimensions[0])-1;i++){ mt[i].Create(dimensions[1]);}
              }
           else{ empty=true; mt=NULL;}
          }

// *******************************************************************

template <class Type,size_t rank>
void TARRAY<Type,rank>::Copy(Type *T)
         { if(Empty()==false)
             { for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ mt[i].Copy(T);} 
              } 
          }

template <class Type,size_t rank>
void TARRAY<Type,rank>::Copy(const TARRAY<Type,rank> &MT)
         { if(Empty()==false)
             { //assert(Dimension(0)==MT.Dimension(0));
               #if !defined(DNDEBUG)
               if(Dimension(0)!=MT.Dimension(0)){ DIFFERENT_DIMENSIONS("Copy(const TARRAY<Type,rank>&","TARRAY<Type,rank>");}
               #endif
               for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ mt[i].Copy(MT.mt[i]);}
              }
          }

// *******************************************************************

// rank >=3
template <class Type,size_t rank> template <size_t rank2>
void TARRAY<Type,rank>::Destroy(TARRAY<Type,rank2>*)
         { if(Empty()==false){ for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ mt[i].Destroy(&mt[i]);}
                               delete []mt;
                               mt=NULL;
                               for(int i=0;i<=static_cast<int>(rank)-1;i++){ dimensions[i]=0;}
                               empty=true;
                              }
          }

// partial specialization for rank 2 tensors
template <class Type,size_t rank> 
void TARRAY<Type,rank>::Destroy(TARRAY<Type,2>*)
         { if(Empty()==false)
             { for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ mt[i].Destroy();} 
               delete []mt;
              }
           mt=NULL;
           empty=true;
          }

// ************************************************************

template <class Type,size_t rank>
TARRAY<Type,rank>::TARRAY(const TARRAY<Type,rank> &MT){ if(MT.Empty()==false){ Create(this,MT.Dimensions()); Copy(MT);} }

template <class Type,size_t rank>
TARRAY<Type,rank>::TARRAY(std::vector<size_t> DIMENSIONS){ if(DIMENSIONS.empty()==false){ Create(this,DIMENSIONS);} }

template <class Type,size_t rank>
TARRAY<Type,rank>::TARRAY(std::vector<size_t> DIMENSIONS,const Type &T){ if(DIMENSIONS.empty()==false){ Create(this,DIMENSIONS); Fill(T);} }

// ************************************************************

template <class Type,size_t rank>
size_t TARRAY<Type,rank>::Dimension(int i) const
         { try{ return dimensions[i];}
           catch(EMPTY &E){ E.Change("Dimension(int)","TARRAY<Type,rank>"); throw E;}
           catch(OUT_OF_RANGE<int> &OOR){ OOR.Change("Dimension(int)","TARRAY<Type,rank>"); throw OOR;}
          }

// *******************************************************************

template <class Type,size_t rank>
void TARRAY<Type,rank>::Fill(const Type &T)
         { for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ (*this)[i].Fill(T);}
          }

//************************************************************************

// Swap index i with index i+1
// rank >=3
template <class Type,size_t rank> template <size_t rank2>
TARRAY<Type,rank>& TARRAY<Type,rank>::SwapAdjacentIndices(TARRAY<Type,rank2>*,int i)
         { //assert(Empty()!=true);
           #if !defined(DNDEBUG)
           if(Empty()==true){ throw EMPTY("SwapAdjacentIndices(TARRAY<Type,rank2>*,int)","TARRAY<Type,rank>");}
           //assert(i>1 && i<rank);
           if(i<=1 || i>=static_cast<int>(rank)){ throw OUT_OF_RANGE<int>(i,1,static_cast<int>(rank),"SwapAdjacentIndices(TARRAY<Type,rank2>*,int)","TARRAY<Type,rank>");}
           #endif

           if(i==1)
             { if(Dimension(0)==Dimension(1))
                 { for(int j=0;j<=static_cast<int>(Dimension(0))-1;j++){ for(int k=j+1;k<=static_cast<int>(Dimension(1))-1;k++){ std::swap((*this)[k][j],(*this)[j][k]);} }
                   return *this;
                  }
               else{ std::vector<size_t> dimensions2(Dimensions());
                     std::swap(dimensions2[0],dimensions2[1]);
                     TARRAY<Type,rank> MT2(dimensions2);
                     for(int j=0;j<=static_cast<int>(Dimension(0))-1;j++){ for(int k=0;k<=static_cast<int>(Dimension(1))-1;k++){ MT2[k][j]=(*this)[j][k];} }
                     return *this=MT2;
                    }
              }
           else{ for(int j=0;j<=static_cast<int>(Dimension(0))-1;j++){ (*this)[j].SwapAdjacentIndices(i-1);}
                 return *this;
                }
          }

// partial specialization for rank 2 tensors
template <class Type,size_t rank> 
TARRAY<Type,2>& TARRAY<Type,rank>::SwapAdjacentIndices(TARRAY<Type,2>*,int i)
         { //assert(Empty()!=true);
           #if !defined(DNDEBUG)
           if(Empty()==true){ throw EMPTY("SwapAdjacentIndices(TARRAY<Type,2>*,int)","TARRAY<Type,2>");}
           //assert(i>=1 && i<=rank);
           if(i<1 || i>rank){ throw OUT_OF_RANGE<int>(i,1,static_cast<int>(rank),"SwapAdjacentIndices(TARRAY<Type,2>*,int)","TARRAY<Type,2>");}
           #endif

           if(i==1)
             { if(Dimension(0)==Dimension(1))
                 { for(int j=0;j<=static_cast<int>(Dimension(0))-1;j++){ for(int k=j+1;k<=static_cast<int>(Dimension(1))-1;k++){ std::swap((*this)[k][j],(*this)[j][k]);} }
                   return *this;
                  }
               else{ std::vector<size_t> dimensions2(Dimensions());
                     std::swap(dimensions2[0],dimensions2[1]);
                     TARRAY<Type,2> MT2(dimensions2);
                     for(int j=0;j<=static_cast<int>(Dimension(0))-1;j++){ for(int k=0;k<=static_cast<int>(Dimension(1))-1;k++){ MT2[k][j]=(*this)[j][k];} }
                     return *this=MT2;
                    }
              }
           else{ return *this;}
          }

// Swap index i with index i+1
template <class Type,size_t rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::SwapAdjacentIndices(int i)
         { try{ return SwapAdjacentIndices(this,i);}
           catch(EMPTY &E){ E.Change("SwapAdjacentIndices(int)","TARRAY<Type,rank>");}
           catch(OUT_OF_RANGE<int> &OOR){ OOR.Change("SwapAdjacentIndices(int)","TARRAY<Type,rank>");}
          }

// *******************************************************************

 // Swap index i with index j
template <class Type,size_t rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::SwapIndices(int i,int j)
         { //assert(Empty()!=true);
           #if !defined(DNDEBUG)
           if(Empty()==true){ throw EMPTY("SwapIndices(int,int)","TARRAY<Type,rank>");}
           //assert(i>=1 && i<=rank); 
           if(i<1 || i>rank){ throw OUT_OF_RANGE<int>(i,1,static_cast<int>(rank),"SwapIndices(int,int)","TARRAY<Type,rank>");}
           //assert(j>=1 && j<=rank);
           if(j<1 || j>rank){ throw OUT_OF_RANGE<int>(j,1,static_cast<int>(rank),"SwapIndices(int,int)","TARRAY<Type,rank>");}
           #endif

           if(i>j){ std::swap(i,j);}

           int i0=i, j0=j;
           while(j>i0+1){ SwapAdjacentIndices(this,j-1); --j;};
           while(i<j0){ SwapAdjacentIndices(this,i); ++i;};

           return *this;
          }

//************************************************************************

template <class Type,size_t rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator=(const TARRAY<Type,rank> &MT)
      { if(MT.Empty()==false){ Destroy(this); Create(this,MT.Dimensions()); Copy(MT);} return (*this);}

template <class Type,size_t rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator+=(const TARRAY<Type,rank> &MT)
         { //assert(Dimension(0)==MT.Dimension(0));
           #if !defined(DNDEBUG)
           if(Dimension(0)!=MT.Dimension(0)){ DIFFERENT_DIMENSIONS("operator+=(const TARRAY<Type,rank>&)","TARRAY<Type,rank>");}
           #endif
           for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ Reference(i)+=MT[i];}
           return (*this);
          }

template <class Type,size_t rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator-=(const TARRAY<Type,rank> &MT)
         { //assert(Dimension(0)==MT.Dimension(0));
           #if !defined(DNDEBUG)
           if(Dimension(0)!=MT.Dimension(0)){ DIFFERENT_DIMENSIONS("operator-=(const TARRAY<Type,rank>&)","TARRAY<Type,rank>");}
           #endif
           for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ Reference(i)-=MT[i];}
           return (*this);
          }

template <class Type,size_t rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator*=(const Type &T)
         { for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ Reference(i)*=T;}
           return (*this);
          }

template <class Type,size_t rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator/=(const Type &T)
         { //assert(T!=Zero<Type>());
           #if !defined(DNDEBUG)
           if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/(const Type&)=","TARRAY<Type,rank>");}
           #endif
           for(int i=0;i<=static_cast<int>(Dimension(0))-1;i++){ Reference(i)/=T;}
           return (*this);
          }

template <class Type,size_t rank>
TARRAY<Type,rank-1>& TARRAY<Type,rank>::operator[](const int &i)
         { //assert(Empty()!=true);
           #if !defined(DNDEBUG)
           if(Empty()==true){ throw EMPTY("TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           //assert(i>=0 && i<=static_cast<int>(Dimension(0))-1);
           if(i<0 || i>static_cast<int>(Dimension(0))-1){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(Dimension(0))-1,"TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           #endif
           return Reference(i);
          }

template <class Type,size_t rank>
TARRAY<Type,rank-1> TARRAY<Type,rank>::operator[](const int &i) const
         { //assert(Empty()!=true);
           #if !defined(DNDEBUG)
           if(Empty()==true){ throw EMPTY("TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           //assert(i>=0 && i<=static_cast<int>(Dimension(0))-1);
           if(i<0 || i>static_cast<int>(Dimension(0))-1){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(Dimension(0))-1,"TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           #endif
           return Value(i);
          }

template <class Type,size_t rank>
Type& TARRAY<Type,rank>::operator[](const std::vector<int> &i0)
         { //assert(static_cast<int>(i0.size())==rank);
           #if !defined(DNDEBUG)
           if(static_cast<int>(i0.size())!=rank){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i0.size()),static_cast<int>(rank),"TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           //assert(i0[0]>=0 && i0[0]<=static_cast<int>(Dimension(0))-1);
           if(i0[0]<0 || i0[0]>static_cast<int>(Dimension(0))-1){ throw OUT_OF_RANGE<int>(i0[0],0,static_cast<int>(Dimension(0))-1,"TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           //assert(Empty()!=true);
           if(Empty()==true){ throw EMPTY("TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           #endif

           std::vector<int> i1(i0.begin()+1,i0.end()); 
           return mt[i0[0]][i1];
           //return Reference(&mt[i0[0]],i0,1);
          }

template <class Type,size_t rank>
Type TARRAY<Type,rank>::operator[](const std::vector<int> &i0) const
         { //assert(static_cast<int>(i0.size())==rank);
           #if !defined(DNDEBUG)
           if(static_cast<int>(i0.size())!=rank){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i0.size()),static_cast<int>(rank),"TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           //assert(i0[0]>=0 && i0[0]<=static_cast<int>(Dimension(0))-1);
           if(i0[0]<0 || i0[0]>static_cast<int>(Dimension(0))-1){ throw OUT_OF_RANGE<int>(i0[0],0,static_cast<int>(Dimension(0))-1,"TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           //assert(Empty()!=true);
           if(Empty()==true){ throw EMPTY("TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           #endif

           std::vector<int> i1(i0.begin()+1,i0.end());
           return mt[i0[0]][i1];
           //return Value(&mt[i0[0]],i0,1);
          }

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

// rank 1 member functions
template <class Type>
void TARRAY<Type,1>::Create(const size_t &DIMENSION)
         { dimension=DIMENSION;

           if(DIMENSION!=0)
             { empty=false;
               mt=new Type[DIMENSION];
              }
           else{ empty=true;
                 mt=NULL;
                }
          }

template <class Type>
void TARRAY<Type,1>::Copy(Type *T)
         { for(int i=0;i<=static_cast<int>(Dimension())-1;i++){ mt[i]=(*T)[i];} 
          }

template <class Type>
void TARRAY<Type,1>::Copy(const std::vector<Type> &v)
         { if(Empty()==false)
             { //assert(Dimension()==v.size());
               #if !defined(DNDEBUG)
               if(Dimension()!=v.size()){ DIFFERENT_DIMENSIONS("Copy(const std::vector<Type>&)","TARRAY<Type,1>");}
               #endif
               for(int i=0;i<=static_cast<int>(Dimension())-1;i++){ mt[i]=v[i];}
              } 
          }

template <class Type>
void TARRAY<Type,1>::Copy(const TARRAY<Type,1> &MT)
         { if(Empty()==false)
             { //assert(Dimension()==MT.Dimension());
               #if !defined(DNDEBUG)
               if(Dimension()!=MT.Dimension()){ DIFFERENT_DIMENSIONS("Copy(const TARRAY<Type,1>&)","TARRAY<Type,1>");}
               #endif
               for(int i=0;i<=static_cast<int>(Dimension())-1;i++){ mt[i]=MT.mt[i];}
              }
          }

template <class Type>
void TARRAY<Type,1>::Destroy(void)
         { if(Empty()==false){ delete []mt;}
           empty=true;
           dimension=0;
           mt=NULL;
          }

//**************************************************************************

template <class Type>
TARRAY<Type,1>::TARRAY(std::vector<int> DIMENSIONS){ if(DIMENSIONS.empty()==false){ Create(DIMENSIONS[0]);} }

template <class Type>
TARRAY<Type,1>::TARRAY(const TARRAY<Type,1> &MT){ if(MT.Empty()==false){ Create(MT.Dimension(0)); Copy(MT);} }

//**************************************************************************

template <class Type>
void TARRAY<Type,1>::Fill(const Type &T)
         { for(int j=0;j<=static_cast<int>(Dimension())-1;j++){ mt[j]=T;}
          }

//**************************************************************************

template <class Type>
size_t TARRAY<Type,1>::Dimension(void) const
         { //assert(Empty()!=true);
           #if !defined(DNDEBUG)
           if(Empty()==true){ throw EMPTY("Dimension(void)","TARRAY<Type,1>");}
           #endif
           return dimension;
          }

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator=(const TARRAY<Type,1> &MT)
      { if(MT.Empty()==false){ Destroy(); Create(MT.Dimension()); Copy(MT);} return (*this);}

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator=(const std::vector<Type> &v)
      { if(v.empty()==false){ Destroy(); Create(v.size()); Copy(v);} return (*this);}

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator+=(const TARRAY<Type,1> &MT)
         { //assert(Dimension()==MT.Dimension());
           #if !defined(DNDEBUG)
           if(Dimension()!=MT.Dimension()){ throw DIFFERENT_DIMENSIONS("operator+=(const TARRAY<Type,1>&)","TARRAY<Type,1>");}
           #endif
           for(int i=0;i<=static_cast<int>(Dimension())-1;i++){ Reference(this,i)+=MT[i];}
           return (*this);
          }

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator-=(const TARRAY<Type,1> &MT)
         { //assert(Dimension()==MT.Dimension());
           #if !defined(DNDEBUG)
           if(Dimension()!=MT.Dimension()){ throw DIFFERENT_DIMENSIONS("operator-=(const TARRAY<Type,1>&)","TARRAY<Type,1>");}
           #endif
           for(int i=0;i<=static_cast<int>(Dimension())-1;i++){ Reference(this,i)-=MT[i];}
           return (*this);
          }

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator*=(const Type &T)
         { for(int i=0;i<=static_cast<int>(Dimension())-1;i++){ Reference(this,i)*=T;}
           return (*this);
          }

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator/=(const Type &T)
         { //assert(T!=Zero<Type>());
           #if !defined(DNDEBUG)
           if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/(const Type&)=","TARRAY<Type,1>");}
           #endif
           for(int i=0;i<=static_cast<int>(Dimension())-1;i++){ Reference(this,i)/=T;}
           return (*this);
          }

template <class Type>
Type& TARRAY<Type,1>::operator[](const int &i)
         { //assert(Empty()!=true);
           #if !defined(DNDEBUG)
           if(Empty()==true){ throw EMPTY("TARRAY<Type,0>& operator[]","TARRAY<Type,1>");}
           if(i<0 || i>static_cast<int>(Dimension())-1){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(Dimension())-1,"TARRAY<Type,0>& operator[]","TARRAY<Type,1>");}
           #endif
           return Reference(i);
          }

template <class Type>
Type TARRAY<Type,1>::operator[](const int &i) const
         { //assert(Empty()!=true);
           #if !defined(DNDEBUG)
           if(Empty()==true){ throw EMPTY("TARRAY<Type,0> operator[]","TARRAY<Type,1>");}
           //assert(i>=0 && i<=static_cast<int>(Dimensions())-1);
           if(i<0 || i>static_cast<int>(Dimension())-1){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(Dimension())-1,"TARRAY<Type,0> operator[]","TARRAY<Type,1>");}
           #endif
           return Value(i);
          }

template <class Type>
Type& TARRAY<Type,1>::operator[](const std::vector<int> &i0)
         { //assert(static_cast<int>(i0.size())==1);
           #if !defined(DNDEBUG)
           if(static_cast<int>(i0.size())!=1){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i0.size()),1,"Type& operator[]","TARRAY<Type,1>");}
           //assert(i0[0]>=0 && i0[0]<=static_cast<int>(Dimensions())-1);
           if(i0[0]<0 || i0[0]>static_cast<int>(Dimension())-1){ throw OUT_OF_RANGE<int>(i0[0],0,static_cast<int>(Dimension())-1,"Type& operator[]","TARRAY<Type,1>");}
           //assert(Empty()!=true); 
           if(Empty()==true){ throw EMPTY("Type& operator[]","TARRAY<Type,1>");}
           #endif
           return Reference(i0[0]);
          }

template <class Type>
Type TARRAY<Type,1>::operator[](const std::vector<int> &i0) const
         { //assert(static_cast<int>(i0.size())==1);
           #if !defined(DNDEBUG)
           if(static_cast<int>(i0.size())!=1){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i0.size()),1,"Type operator[]","TARRAY<Type,1>");}
           //assert(i0[0]>=0 && i0[0]<=static_cast<int>(Dimensions())-1);
           if(i0[0]<0 || i0[0]>static_cast<int>(Dimension())-1){ throw OUT_OF_RANGE<int>(i0[0],0,static_cast<int>(Dimension())-1,"Type operator[]","TARRAY<Type,1>");}
           //assert(Empty()!=true);
           if(Empty()==true){ throw EMPTY("Type operator[]","TARRAY<Type,1>");}
           #endif
           return Value(i0[0]);
          }

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

#endif
