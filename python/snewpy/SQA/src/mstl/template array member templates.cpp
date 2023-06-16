#if !defined(_TARRAY_MEMBERS)
#define _TARRAY_MEMBERS

// *******************************************************************
// *******************************************************************
// *******************************************************************

// rank >=3
template <class Type,int rank> template <int rank2>
void TARRAY<Type,rank>::Create(TARRAY<Type,rank2>*,std::vector<int> DIMENSIONS)
         { if(DIMENSIONS.size()!=rank){ throw DIFFERENT_DIMENSIONS("Create(TARRAY<Type,rank2>*,std::vector<int>)","TARRAY<Type,rank>");}
           dimensions=DIMENSIONS;

           if(dimensions[0]!=0)
             { empty=false;
               mt=new TARRAY<Type,rank-1>[dimensions[0]];
               DIMENSIONS.erase(DIMENSIONS.begin());

               int i;
               for(i=0;i<=dimensions[0]-1;i++){ mt[i].Create(&mt[i],DIMENSIONS);}
              }
           else{ empty=true; mt=NULL;}
          }

// partial specialization for rank 2 tensors
template <class Type,int rank> 
void TARRAY<Type,rank>::Create(TARRAY<Type,2>*,std::vector<int> DIMENSIONS)
         { if(DIMENSIONS.size()!=2){ throw DIFFERENT_DIMENSIONS("Create(TARRAY<Type,2>*,std::vector<int>)","TARRAY<Type,1>");}
           dimensions=DIMENSIONS;

           if(dimensions[0]!=0)
             { empty=false;
               mt=new TARRAY<Type,1>[dimensions[0]];

               int i;
               for(i=0;i<=dimensions[0]-1;i++){ mt[i].Create(dimensions[1]);}
              }
           else{ empty=true; mt=NULL;}
          }

// *******************************************************************

template <class Type,int rank>
void TARRAY<Type,rank>::Copy(Type *T)
         { if(Empty()==false)
             { int i;
               for(i=0;i<=Dimension(0)-1;i++){ mt[i].Copy(T);} 
              } 
          }

template <class Type,int rank>
void TARRAY<Type,rank>::Copy(const TARRAY<Type,rank> &MT)
         { if(Empty()==false)
             { if(Dimension(0)!=MT.Dimension(0)){ DIFFERENT_DIMENSIONS("Copy(const TARRAY<Type,rank>&","TARRAY<Type,rank>");}
               int i;
               for(i=0;i<=Dimension(0)-1;i++){ mt[i].Copy(MT.mt[i]);}
              }
          }

// *******************************************************************

// rank >=3
template <class Type,int rank> template <int rank2>
void TARRAY<Type,rank>::Destroy(TARRAY<Type,rank2>*)
         { if(Empty()==false){ int i;
                               for(i=0;i<=Dimension(0)-1;i++){ mt[i].Destroy(&mt[i]);}
                               delete []mt;

                               mt=NULL;

                               for(i=0;i<=rank-1;i++){ dimensions[i]=0;}
                               empty=true;
                              }
          }

// partial specialization for rank 2 tensors
template <class Type,int rank> 
void TARRAY<Type,rank>::Destroy(TARRAY<Type,2>*)
         { if(Empty()==false)
             { int i;
               for(i=0;i<=Dimension(0)-1;i++){ mt[i].Destroy();} 
               delete []mt;
              }
           mt=NULL;
           empty=true;
          }

// ************************************************************

template <class Type,int rank>
TARRAY<Type,rank>::TARRAY(const TARRAY<Type,rank> &MT){ if(MT.Empty()==false){ Create(this,MT.Dimensions()); Copy(MT);} }

template <class Type,int rank>
TARRAY<Type,rank>::TARRAY(std::vector<int> DIMENSIONS){ if(DIMENSIONS.empty()==false){ Create(this,DIMENSIONS);} }

template <class Type,int rank>
TARRAY<Type,rank>::TARRAY(std::vector<int> DIMENSIONS,const Type &T){ if(DIMENSIONS.empty()==false){ Create(this,DIMENSIONS); Fill(T);} }

// ************************************************************

template <class Type,int rank>
int TARRAY<Type,rank>::Dimension(int i) const
         { try{ return dimensions[i];}
           catch(EMPTY &E){ E.Change("Dimension(int)","TARRAY<Type,rank>"); throw E;}
           catch(OUT_OF_RANGE<int> &OOR){ OOR.Change("Dimension(int)","TARRAY<Type,rank>"); throw OOR;}
          }

// *******************************************************************

template <class Type,int rank>
void TARRAY<Type,rank>::Fill(const Type &T)
         { for(int i=0;i<=Dimension(0)-1;i++){ (*this)[i].Fill(T);}
          }

//************************************************************************

// Swap index i with index i+1
// rank >=3
template <class Type,int rank> template <int rank2>
TARRAY<Type,rank>& TARRAY<Type,rank>::SwapAdjacentIndices(TARRAY<Type,rank2>*,int i)
         { if(Empty()==true){ throw EMPTY("SwapAdjacentIndices(TARRAY<Type,rank2>*,int)","TARRAY<Type,rank>");}
           if(i<1 || i>rank){ throw OUT_OF_RANGE<int>(i,1,rank,"SwapAdjacentIndices(TARRAY<Type,rank2>*,int)","TARRAY<Type,rank>");}

           if(i==1)
             { if(Dimension(0)==Dimension(1))
                 { int j,k;
                   for(j=0;j<=Dimension(0)-1;j++){ for(k=j+1;k<=Dimension(1)-1;k++){ std::swap((*this)[k][j],(*this)[j][k]);} }
                   return *this;
                  }
               else{ std::vector<int> dimensions2(Dimensions());
                     std::swap(dimensions2[0],dimensions2[1]);
                     TARRAY<Type,rank> MT2(dimensions2);
                     int j,k;
                     for(j=0;j<=Dimension(0)-1;j++){ for(k=0;k<=Dimension(1)-1;k++){ MT2[k][j]=(*this)[j][k];} }
                     return *this=MT2;
                    }
              }
           else{ int j;
                 for(j=0;j<=Dimension(0)-1;j++){ (*this)[j].SwapAdjacentIndices(i-1);}
                 return *this;
                }
          }

// partial specialization for rank 2 tensors
template <class Type,int rank> 
TARRAY<Type,2>& TARRAY<Type,rank>::SwapAdjacentIndices(TARRAY<Type,2>*,int i)
         { if(Empty()==true){ throw EMPTY("SwapAdjacentIndices(TARRAY<Type,2>*,int)","TARRAY<Type,2>");}
           if(i<1 || i>rank){ throw OUT_OF_RANGE<int>(i,1,rank,"SwapAdjacentIndices(TARRAY<Type,2>*,int)","TARRAY<Type,2>");}

           if(i==1)
             { if(Dimension(0)==Dimension(1))
                 { int j,k;
                   for(j=0;j<=Dimension(0)-1;j++){ for(k=j+1;k<=Dimension(1)-1;k++){ std::swap((*this)[k][j],(*this)[j][k]);} }
                   return *this;
                  }
               else{ std::vector<int> dimensions2(Dimensions());
                     std::swap(dimensions2[0],dimensions2[1]);
                     TARRAY<Type,2> MT2(dimensions2);
                     int j,k;
                     for(j=0;j<=Dimension(0)-1;j++){ for(int k=0;k<=Dimension(1)-1;k++){ MT2[k][j]=(*this)[j][k];} }
                     return *this=MT2;
                    }
              }
           else{ return *this;}
          }

// Swap index i with index i+1
template <class Type,int rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::SwapAdjacentIndices(int i)
         { try{ return SwapAdjacentIndices(this,i);}
           catch(EMPTY &E){ E.Change("SwapAdjacentIndices(int)","TARRAY<Type,rank>");}
           catch(OUT_OF_RANGE<int> &OOR){ OOR.Change("SwapAdjacentIndices(int)","TARRAY<Type,rank>");}
          }

// *******************************************************************

 // Swap index i with index j
template <class Type,int rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::SwapIndices(int i,int j)
         { if(Empty()==true){ throw EMPTY("SwapIndices(int,int)","TARRAY<Type,rank>");}
           if(i>j){ std::swap(i,j);}
           if(i<1 || i>rank){ throw OUT_OF_RANGE<int>(i,1,rank,"SwapIndices(int,int)","TARRAY<Type,rank>");}
           if(j<1 || j>rank){ throw OUT_OF_RANGE<int>(j,1,rank,"SwapIndices(int,int)","TARRAY<Type,rank>");}

           int i0=i, j0=j;
           while(j>i0+1){ SwapAdjacentIndices(this,j-1); --j;};
           while(i<j0){ SwapAdjacentIndices(this,i); ++i;};

           return *this;
          }

//************************************************************************

template <class Type,int rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator=(const TARRAY<Type,rank> &MT)
      { if(MT.Empty()==false){ Destroy(this); Create(this,MT.Dimensions()); Copy(MT);} return (*this);}

template <class Type,int rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator+=(const TARRAY<Type,rank> &MT)
         { if(Dimension(0)!=MT.Dimension(0)){ DIFFERENT_DIMENSIONS("operator+=(const TARRAY<Type,rank>&)","TARRAY<Type,rank>");}
           int i;
           for(i=0;i<=Dimension(0)-1;i++){ Reference(i)+=MT[i];}
           return (*this);
          }

template <class Type,int rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator-=(const TARRAY<Type,rank> &MT)
         { if(Dimension(0)!=MT.Dimension(0)){ DIFFERENT_DIMENSIONS("operator-=(const TARRAY<Type,rank>&)","TARRAY<Type,rank>");}
           int i;
           for(i=0;i<=Dimension(0)-1;i++){ Reference(i)-=MT[i];}
           return (*this);
          }

template <class Type,int rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator*=(const Type &T)
         { int i;
           for(i=0;i<=Dimension(0)-1;i++){ Reference(i)*=T;}
           return (*this);
          }

template <class Type,int rank>
TARRAY<Type,rank>& TARRAY<Type,rank>::operator/=(const Type &T)
         { if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/(const Type&)=","TARRAY<Type,rank>");}
           int i;
           for(i=0;i<=Dimension(0)-1;i++){ Reference(i)/=T;}
           return (*this);
          }

template <class Type,int rank>
TARRAY<Type,rank-1>& TARRAY<Type,rank>::operator[](const int &i)
         { if(Empty()==true){ throw EMPTY("TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           if(i<0 || i>Dimension(0)-1){ throw OUT_OF_RANGE<int>(i,0,Dimension(0)-1,"TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           return Reference(i);
          }

template <class Type,int rank>
TARRAY<Type,rank-1> TARRAY<Type,rank>::operator[](const int &i) const
         { if(Empty()==true){ throw EMPTY("TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           if(i<0 || i>Dimension(0)-1){ throw OUT_OF_RANGE<int>(i,0,Dimension(0)-1,"TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           return Value(i);
          }

template <class Type,int rank>
Type& TARRAY<Type,rank>::operator[](const std::vector<int> &i0)
         { if(static_cast<int>(i0.size())!=rank){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i0.size()),rank,"TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           if(i0[0]<0 || i0[0]>Dimension(0)-1){ throw OUT_OF_RANGE<int>(i0[0],0,Dimension(0)-1,"TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           if(Empty()==true){ throw EMPTY("TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           std::vector<int> i1(i0.begin()+1,i0.end()); 
           return mt[i0[0]][i1];
           //return Reference(&mt[i0[0]],i0,1);
          }

template <class Type,int rank>
Type TARRAY<Type,rank>::operator[](const std::vector<int> &i0) const
         { if(static_cast<int>(i0.size())!=rank){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i0.size()),rank,"TARRAY<Type,rank-1>& operator[]","TARRAY<Type,rank>");}
           if(i0[0]<0 || i0[0]>Dimension(0)-1){ throw OUT_OF_RANGE<int>(i0[0],0,Dimension(0)-1,"TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           if(Empty()==true){ throw EMPTY("TARRAY<Type,rank-1> operator[]","TARRAY<Type,rank>");}
           std::vector<int> i1(i0.begin()+1,i0.end());
           return mt[i0[0]][i1];
           //return Value(&mt[i0[0]],i0,1);
          }

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

// rank 1 member functions
template <class Type>
void TARRAY<Type,1>::Create(const int &DIMENSION)
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
         { int i;
           for(i=0;i<=Dimension()-1;i++){ mt[i]=(*T)[i];} 
          }

template <class Type>
void TARRAY<Type,1>::Copy(const std::vector<Type> &v)
         { if(Empty()==false)
             { if(Dimension()!=v.size()){ DIFFERENT_DIMENSIONS("Copy(const std::vector<Type>&)","TARRAY<Type,1>");}
               int i;
               for(i=0;i<=Dimension()-1;i++){ mt[i]=v[i];}
              } 
          }

template <class Type>
void TARRAY<Type,1>::Copy(const TARRAY<Type,1> &MT)
         { if(Empty()==false)
             { if(Dimension()!=MT.Dimension()){ DIFFERENT_DIMENSIONS("Copy(const TARRAY<Type,1>&)","TARRAY<Type,1>");}
               int i;
               for(i=0;i<=Dimension()-1;i++){ mt[i]=MT.mt[i];}
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
         { for(int j=0;j<=Dimension()-1;j++){ mt[j]=T;}
          }

//**************************************************************************

template <class Type>
int TARRAY<Type,1>::Dimension(void) const
         { if(Empty()==true){ throw EMPTY("Dimension(void)","TARRAY<Type,1>");}
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
         { if(Dimension()!=MT.Dimension()){ throw DIFFERENT_DIMENSIONS("operator+=(const TARRAY<Type,1>&)","TARRAY<Type,1>");}
           int i;
           for(i=0;i<=Dimension()-1;i++){ Reference(this,i)+=MT[i];}
           return (*this);
          }

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator-=(const TARRAY<Type,1> &MT)
         { if(Dimension()!=MT.Dimension()){ throw DIFFERENT_DIMENSIONS("operator-=(const TARRAY<Type,1>&)","TARRAY<Type,1>");}
           int i;
           for(i=0;i<=Dimension()-1;i++){ Reference(this,i)-=MT[i];}
           return (*this);
          }

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator*=(const Type &T)
         { int i;
           for(i=0;i<=Dimension()-1;i++){ Reference(this,i)*=T;}
           return (*this);
          }

template <class Type>
TARRAY<Type,1>& TARRAY<Type,1>::operator/=(const Type &T)
         { if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/(const Type&)=","TARRAY<Type,1>");}
           int i;
           for(i=0;i<=Dimension()-1;i++){ Reference(this,i)/=T;}
           return (*this);
          }

template <class Type>
Type& TARRAY<Type,1>::operator[](const int &i)
         { if(Empty()==true){ throw EMPTY("TARRAY<Type,0>& operator[]","TARRAY<Type,1>");}
           if(i<0 || i>Dimension()-1){ throw OUT_OF_RANGE<int>(i,0,Dimension()-1,"TARRAY<Type,0>& operator[]","TARRAY<Type,1>");}
           return Reference(i);
          }

template <class Type>
Type TARRAY<Type,1>::operator[](const int &i) const
         { if(Empty()==true){ throw EMPTY("TARRAY<Type,0> operator[]","TARRAY<Type,1>");}
           if(i<0 || i>Dimension()-1){ throw OUT_OF_RANGE<int>(i,0,Dimension()-1,"TARRAY<Type,0> operator[]","TARRAY<Type,1>");}
           return Value(i);
          }

template <class Type>
Type& TARRAY<Type,1>::operator[](const std::vector<int> &i0)
         { if(static_cast<int>(i0.size())!=1){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i0.size()),1,"Type& operator[]","TARRAY<Type,1>");}
           if(i0[0]<0 || i0[0]>Dimension()-1){ throw OUT_OF_RANGE<int>(i0[0],0,Dimension()-1,"Type& operator[]","TARRAY<Type,1>");}
           if(Empty()==true){ throw EMPTY("Type& operator[]","TARRAY<Type,1>");}
           return Reference(i0[0]);
          }

template <class Type>
Type TARRAY<Type,1>::operator[](const std::vector<int> &i0) const
         { if(static_cast<int>(i0.size())!=1){ throw NOT_EQUAL_VALUES<int>(static_cast<int>(i0.size()),1,"Type operator[]","TARRAY<Type,1>");}
           if(i0[0]<0 || i0[0]>Dimension()-1){ throw OUT_OF_RANGE<int>(i0[0],0,Dimension()-1,"Type operator[]","TARRAY<Type,1>");}
           if(Empty()==true){ throw EMPTY("Type operator[]","TARRAY<Type,1>");}
           return Value(i0[0]);
          }

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

#endif
