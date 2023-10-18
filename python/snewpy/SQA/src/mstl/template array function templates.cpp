#if !defined(_TARRAY_FUNCTIONS)
#define _TARRAY_FUNCTIONS

template <class Type,size_t rank1,size_t rank2>
TARRAY<Type,rank1+rank2> operator*(const TARRAY<Type,rank1> &MT1,const TARRAY<Type,rank2> &MT2)
         { TARRAY<Type,rank1+rank2> MT3(MT1.Dimensions()+MT2.Dimensions());
           int i;
           for(i=0;i<=static_cast<int>(MT1.Dimension())-1;i++){ MT3[i]=MT1[i]*MT2;}
           return MT3;
          }

// partial specialization for zeroth rank tensors
template <class Type,size_t rank1>
TARRAY<Type,rank1> operator*(const TARRAY<Type,rank1> &MT1,const TARRAY<Type,0> &MT2)
         { return MT1*(Type)MT2;}

template <class Type,size_t rank2>
TARRAY<Type,rank2> operator*(const TARRAY<Type,0> &MT1,const TARRAY<Type,rank2> &MT2)
         { return (Type)MT1*MT2;} 

template <class Type>
TARRAY<Type,0> operator*(const TARRAY<Type,0> &MT1,const TARRAY<Type,0> &MT2)
         { return static_cast<Type>(MT1)*static_cast<Type>(MT2);} 
 
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

template <class Type,size_t rank>
TARRAY<Type,rank> operator*(const TARRAY<Type,rank> &MT1,const Type &T)
         { TARRAY<Type,rank> MT2(MT1.Dimensions());
           for(int i=0;i<=static_cast<int>(MT2.Dimension())-1;i++){ MT2[i]=MT1[i]*T;}
           return MT2;
          }

// partial specialization for first rank tensors
template <class Type>
TARRAY<Type,1> operator*(const TARRAY<Type,1> &MT1,const Type &T)
         { TARRAY<Type,1> MT2(MT1.Dimension());
           for(int i=0;i<=static_cast<int>(MT2.Dimension())-1;i++){ MT2[i]=MT1[i]*T;}
           return MT2;
          }

// partial specialization for zeroth rank tensors
template <class Type>
TARRAY<Type,0> operator*(const TARRAY<Type,0> &MT,const Type &T)
         { return static_cast<Type>(MT)*T;}
 
//*********************************************************************************** 
//***********************************************************************************
//*********************************************************************************** 

template <class Type,size_t rank>
TARRAY<Type,rank> operator*(const Type &T,const TARRAY<Type,rank> &MT1)
         { TARRAY<Type,rank> MT2(MT1.Dimensions()); 
           for(int i=0;i<=static_cast<int>(MT1.Dimension())-1;i++){ MT2[i]=T*MT1[i];} 
           return MT2;
          }

// partial specialization for first rank tensors
template <class Type>
TARRAY<Type,1> operator*(const Type &T,const TARRAY<Type,1> &MT1)
         { TARRAY<Type,1> MT2(MT1.Dimension());
           for(int i=0;i<=static_cast<int>(MT1.Dimension())-1;i++){ MT2[i]=T*MT1[i];}
           return MT2;
          }

// partial specialization for zeroth rank tensors
template <class Type>
TARRAY<Type,0> operator*(const Type &T,const TARRAY<Type,0> &MT)
         { return T*static_cast<Type>(MT);}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

template <class Type,size_t rank>
TARRAY<Type,rank> operator+(const TARRAY<Type,rank> &MT1,const TARRAY<Type,rank> &MT2)
         { TARRAY<Type,rank> MT3(MT1.Dimensions());
           for(int i=0;i<=static_cast<int>(MT1.Dimension())-1;i++){ MT3[i]=MT1[i]+MT2[i];}
           return MT3;
          }

// partial specialization for zeroth rank tensors
template <class Type>
TARRAY<Type,0> operator+(const TARRAY<Type,0> &MT1,const TARRAY<Type,0> &MT2)
         { return static_cast<Type>(MT1)+static_cast<Type>(MT2);}

//***********************************************************************************
//***********************************************************************************
//*********************************************************************************** 

template <class Type,size_t rank>
TARRAY<Type,rank> operator-(const TARRAY<Type,rank> &MT1,const TARRAY<Type,rank> &MT2)
         { TARRAY<Type,rank> MT3(MT1.Dimensions()); 
           for(int i=0;i<=static_cast<int>(MT1.Dimension())-1;i++){ MT3[i]=MT1[i]-MT2[i];}
           return MT3; 
          } 
 
// partial specialization for zeroth rank tensors 
template <class Type>
TARRAY<Type,0> operator-(const TARRAY<Type,0> &MT1,const TARRAY<Type,0> &MT2)
         { return static_cast<Type>(MT1)-static_cast<Type>(MT2);} 

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

// swaps index i with index i+1
template <class Type,size_t rank>
TARRAY<Type,rank> SwapAdjacentIndices(const TARRAY<Type,rank>& MT,int i)
         { if(i==1){ std::vector<size_t> dimensions2(MT.Dimensions());
                     std::swap(dimensions2[0],dimensions2[1]);
                     TARRAY<Type,rank> MT2(dimensions2);
                     for(int j=0;j<=static_cast<int>(MT.Dimension(0))-1;j++){ for(int k=0;k<=static_cast<int>(MT.Dimension(1))-1;k++){ MT2[k][j]=MT[j][k];} }
                     return MT2;
                    }
           else{ TARRAY<Type,rank> MT2(MT);
                 for(int j=0;j<=static_cast<int>(MT.Dimension(0))-1;j++){ MT2[j]=SwapAdjacentIndices(MT2[j],i-1);}
                 return MT2;
                }
          }

template <class Type>
TARRAY<Type,2> SwapAdjacentIndices(const TARRAY<Type,2>& MT,int)
         { std::vector<size_t> dimensions2(MT.Dimensions());
           std::swap(dimensions2[0],dimensions2[1]);
           TARRAY<Type,2> MT2(dimensions2);
           for(int j=0;j<=static_cast<int>(MT.Dimension(0))-1;j++){ for(int k=0;k<=static_cast<int>(MT.Dimension(1))-1;k++){ MT2[k][j]=MT[j][k];} }
           return MT2;
          }

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

template <class Type,size_t rank>
TARRAY<Type,rank-2> Contract(const TARRAY<Type,rank> &MT,int i,int j)
         { if(i>j){ std::swap(i,j);}
           if(MT.Dimension(i-1)!=MT.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,rank>&,int,int)");}

           if(i==1 && j==2)
             { std::vector<size_t> dimensions1(MT.Dimensions());
               dimensions1.erase(dimensions1.begin());
               dimensions1.erase(dimensions1.begin());

               TARRAY<Type,rank-2> MT1(dimensions1);

               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ MT1[k]=MT[0][0][k];}

               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++)
                  { for(int l=1;l<=static_cast<int>(MT.Dimension(0))-1;l++){ MT1[k]+=MT[l][l][k];} 
                   }
               return MT1;
             }
           else{ TARRAY<Type,rank> MT1(MT);
                 MT1.SwapIndices(1,i);
                 MT1.SwapIndices(2,j);
                 return Contract(MT1,1,2);
                }
          }

// partial specialization for rank 2 tensors
template <class Type>
TARRAY<Type,0> Contract(const TARRAY<Type,2> &MT,int,int)
         { if(MT.Dimension(0)!=MT.Dimension(1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,2>&,int,int)");}
           TARRAY<Type,0> T(MT[0][0]);
           for(int i=1;i<=static_cast<int>(MT.Dimension())-1;i++){ T+=MT[i][i];}
           return T;
          }

// partial specialization for rank 3 tensors
template <class Type>
TARRAY<Type,1> Contract(const TARRAY<Type,3> &MT,int i,int j)
         { //if(i>j){ std::swap(i,j);}
           if(MT.Dimension(i-1)!=MT.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,rank>&,int,int)");}

           if( (i==1 && j==2) || (i==2 && j==1) )
             { TARRAY<Type,1> MT1(MT.Dimension(2));
               for(int k=0;k<=static_cast<int>(MT1.Dimension())-1;k++){ MT1[k]=MT[0][0][k];}
               for(int k=0;k<=static_cast<int>(MT1.Dimension())-1;k++){ for(int l=1;l<=static_cast<int>(MT.Dimension(0))-1;l++){ MT1[k]+=MT[l][l][k];} }
               return MT1;
              }
           if( (i==1 && j==3) || (i==2 && j==1) )
             { TARRAY<Type,1> MT1(MT.Dimension(1));
               for(int k=0;k<=static_cast<int>(MT1.Dimension())-1;k++){ MT1[k]=MT[0][k][0];}
               for(int k=0;k<=static_cast<int>(MT1.Dimension())-1;k++){ for(int l=1;l<=MT.Dimension(2)-1;l++){ MT1[k]+=MT[l][k][l];} }
               return MT1;
              }
           if( (i==2 && j==3) || (i==3 && j==2) )
             { TARRAY<Type,1> MT1(MT.Dimension(0));
               for(int k=0;k<=static_cast<int>(MT1.Dimension())-1;k++){ MT1[k]=MT[k][0][0];}
               for(int k=0;k<=static_cast<int>(MT1.Dimension())-1;k++){ for(int l=1;l<=static_cast<int>(MT.Dimension(1))-1;l++){ MT1[k]+=MT[k][l][l];} }
               return MT1;
              }
          }

// partial specialization for rank 4 tensors
template <class Type>
TARRAY<Type,2> Contract(const TARRAY<Type,4> &MT,int i,int j)
         { //if(i>j){ std::swap(i,j);}
           if(MT.Dimension(i-1)!=MT.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,rank>&,int,int)");}

           std::vector<size_t> dimensions(2); 
           if( (i==1 && j==2) || (i==2 && j==1) )
             { dimensions[0]=MT.Dimension(2); dimensions[1]=MT.Dimension(3);
               TARRAY<Type,2> MT1(dimensions);
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ MT1[k][l]=MT[0][0][k][l];} }
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT.Dimension(0))-1;m++){ MT1[k][l]+=MT[m][m][k][l];} } }
               return MT1;
              }
           if( (i==1 && j==3) || (i==3 && j==1) )
             { dimensions[0]=MT.Dimension(1); dimensions[1]=MT.Dimension(3);
               TARRAY<Type,2> MT1(dimensions);
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ MT1[k][l]=MT[0][k][0][l];} }
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT.Dimension(0))-1;m++){ MT1[k][l]+=MT[m][k][m][l];} } }
               return MT1;
              }
           if( (i==1 && j==4) || (i==4 && j==1) )
             { dimensions[0]=MT.Dimension(1); dimensions[1]=MT.Dimension(2);
               TARRAY<Type,2> MT1(dimensions);
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ MT1[k][l]=MT[0][k][l][0];} }
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT.Dimension(0))-1;m++){ MT1[k][l]+=MT[m][k][l][m];} } }
               return MT1;
              }
           if( (i==2 && j==3) || (i==3 && j==2) )
             { dimensions[0]=MT.Dimension(0); dimensions[1]=MT.Dimension(3);
               TARRAY<Type,2> MT1(dimensions);
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ MT1[k][l]=MT[k][0][0][l];} }
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT.Dimension(0))-1;m++){ MT1[k][l]+=MT[k][m][m][l];} } }
               return MT1;
              }
           if( (i==2 && j==4) || (i==4 && j==2) )
             { dimensions[0]=MT.Dimension(0); dimensions[1]=MT.Dimension(2);
               TARRAY<Type,2> MT1(dimensions);
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ MT1[k][l]=MT[k][0][l][0];} }
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT.Dimension(0))-1;m++){ MT1[k][l]+=MT[k][m][l][m];} } }
               return MT1;
              }
           if( (i==3 && j==4) || (i==4 && j==3) )
             { dimensions[0]=MT.Dimension(0); dimensions[1]=MT.Dimension(1);
               TARRAY<Type,2> MT1(dimensions);
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ MT1[k][l]=MT[k][l][0][0];} }
               for(int k=0;k<=static_cast<int>(MT1.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT1.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT.Dimension(0))-1;m++){ MT1[k][l]+=MT[k][l][m][m];} } }
               return MT1;
              }
          }

//***********************************************************************************

// This is the general function. There are many specializations for this function so that low rank tensors, 0,1,2,3,4
// are dealt with more quickly than in the general case.
template <class Type,size_t rank1,size_t rank2>
TARRAY<Type,rank1+rank2-2> Contract(const TARRAY<Type,rank1> &MT1,const TARRAY<Type,rank2> &MT2,int i,int j)
         { if(MT1.Dimension(i-1)!=MT2.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type>&,const TARRAY<Type>&,int,int)");}

           if(i==1 && j==1)
             { std::vector<size_t> dimensions3(MT1.Dimensions());
               dimensions3.erase(dimensions3.begin());
               for(int k=0;k<=static_cast<int>(MT2.Rank())-1;k++){ dimensions3.push_back(MT2.Dimension(k));}
               dimensions3.erase(dimensions3.begin()+rank1);

               TARRAY<Type,rank1+rank2-2> MT3(dimensions3);

               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ MT3[k]=MT1[0][k]*MT2[0];}
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=1;l<=static_cast<int>(MT1.Dimension(0))-1;l++){ MT3[k]+=MT1[l][k]*MT2[l];} }
               return MT3;
              }
           else{ TARRAY<Type,rank1> MT3(MT1); TARRAY<Type,rank2> MT4(MT2);
                 MT3.SwapIndices(i,1);
                 MT4.SwapIndices(j,1);
                 return Contract(MT3,MT4,1,1);
                }
          }

//partial specialization
template <class Type,size_t rank2>
TARRAY<Type,1+rank2-2> Contract(const TARRAY<Type,1> &MT1,const TARRAY<Type,rank2> &MT2,int,int j)
         { try{ return Contract(MT1,MT2,j);}
           catch(DIFFERENT_DIMENSIONS &DD){ DD.Change("Contract(const TARRAY<Type,1>&,const TARRAY<Type>&,int,int)"); throw DD;}
          }

//partial specialization
template <class Type,size_t rank1>
TARRAY<Type,rank1+1-2> Contract(const TARRAY<Type,rank1> &MT1,const TARRAY<Type,1> &MT2,int i,int)
         { try{ return Contract(MT1,MT2,i);}
           catch(DIFFERENT_DIMENSIONS &DD){ DD.Change("Contract(const TARRAY<Type>&,const TARRAY<Type,1>&,int,int)"); throw DD;}
          }

//partial specialization
template <class Type>
TARRAY<Type,0> Contract(const TARRAY<Type,1> &MT1,const TARRAY<Type,1> &MT2,int,int)
         { try{ return Contract(MT1,MT2);}
           catch(DIFFERENT_DIMENSIONS &DD){ DD.Change("Contract(const TARRAY<Type,1>&,const TARRAY<Type,1>&,int,int)"); throw DD;}
          }

//partial specialization for second rank tensors
template <class Type>
TARRAY<Type,2> Contract(const TARRAY<Type,2> &MT1,const TARRAY<Type,2> &MT2,int i,int j)
         { if(MT1.Dimension(i-1)!=MT2.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,2>&,const TARRAY<Type,2>&,int,int)");}

           std::vector<size_t> dimensions3(MT1.Dimensions());
           dimensions3.erase(dimensions3.begin()+i-1);
           for(int k=0;k<=static_cast<int>(MT2.Rank())-1;k++){ dimensions3.push_back(MT2.Dimension(k));}
           dimensions3.erase(dimensions3.begin()+j);   // j=(2-1) + (j-1)

           TARRAY<Type,2> MT3(dimensions3);

           if(i==1 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ MT3[k][l]=MT1[0][k]*MT2[0][l];} }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT1.Dimension(0))-1;m++){ MT3[k][l]+=MT1[m][k]*MT2[m][l];} } }
              }
           if(i==1 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ MT3[k][l]=MT1[0][k]*MT2[l][0];} }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT1.Dimension(0))-1;m++){ MT3[k][l]+=MT1[m][k]*MT2[l][m];} } }
              }
           if(i==2 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ MT3[k][l]=MT1[k][0]*MT2[0][l];} }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT1.Dimension(1))-1;m++){ MT3[k][l]+=MT1[k][m]*MT2[m][l];} } }
              }
           if(i==2 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ MT3[k][l]=MT1[k][0]*MT2[l][0];} }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=1;m<=static_cast<int>(MT1.Dimension(1))-1;m++){ MT3[k][l]+=MT1[k][m]*MT2[l][m];} } }
              }

           return MT3;
          }

//partial specialization for second rank tensor and third rank
template <class Type>
TARRAY<Type,3> Contract(const TARRAY<Type,2> &MT1,const TARRAY<Type,3> &MT2,int i,int j)
         { if(MT1.Dimension(i-1)!=MT2.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,2>&,const TARRAY<Type,3>&,int,int)");}

           std::vector<int> dimensions3(MT1.Dimensions());
           dimensions3.erase(dimensions3.begin()+i-1);
           for(int k=0;k<=static_cast<int>(MT2.Rank())-1;k++){ dimensions3.push_back(MT2.Dimension(k));}
           dimensions3.erase(dimensions3.begin()+j);    // j=(2-1) + (j-1)

           TARRAY<Type,3> MT3(dimensions3);

           if(i==1 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[0][k]*MT2[0][l][m];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(0))-1;n++){ MT3[k][l][m]+=MT1[n][k]*MT2[n][l][m];} } } }
              }
           if(i==1 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[0][k]*MT2[l][0][m];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(0))-1;n++){ MT3[k][l][m]+=MT1[n][k]*MT2[l][n][m];} } } }
              }
           if(i==1 && j==3)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[0][k]*MT2[l][m][0];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(0))-1;n++){ MT3[k][l][m]+=MT1[n][k]*MT2[l][m][n];} } } }
              }

           if(i==2 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[k][0]*MT2[0][l][m];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(1))-1;n++){ MT3[k][l][m]+=MT1[k][n]*MT2[n][l][m];} } } }
              }
           if(i==2 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[k][0]*MT2[l][0][m];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(1))-1;n++){ MT3[k][l][m]+=MT1[k][n]*MT2[l][n][m];} } } }
              }
           if(i==2 && j==3)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[k][0]*MT2[l][m][0];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(1))-1;n++){ MT3[k][l][m]+=MT1[k][n]*MT2[l][m][n];} } } }
              }

           return MT3;
          }

//partial specialization for third rank tensor and second rank
template <class Type>
TARRAY<Type,3> Contract(const TARRAY<Type,3> &MT1,const TARRAY<Type,2> &MT2,int i,int j)
         { if(MT1.Dimension(i-1)!=MT2.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,2>&,const TARRAY<Type,3>&,int,int)");}

           std::vector<size_t> dimensions3(MT1.Dimensions());
           dimensions3.erase(dimensions3.begin()+i-1);

           for(int k=0;k<=static_cast<int>(MT2.Rank())-1;k++){ dimensions3.push_back(MT2.Dimension(k));}
           dimensions3.erase(dimensions3.begin()+1+j);  // 1+j=(3-1) + (j-1)

           TARRAY<Type,3> MT3(dimensions3);

           if(i==1 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[0][k][l]*MT2[0][m];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(0))-1;n++){ MT3[k][l][m]+=MT1[n][k][l]*MT2[n][m];} } } }
              }
           if(i==1 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[0][k][l]*MT2[m][0];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(0))-1;n++){ MT3[k][l][m]+=MT1[n][k][l]*MT2[m][n];} } } }
              }

           if(i==2 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[k][0][l]*MT2[0][m];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(1))-1;n++){ MT3[k][l][m]+=MT1[k][n][l]*MT2[n][m];} } } }
              }
           if(i==2 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[k][0][l]*MT2[m][0];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=static_cast<int>(MT1.Dimension(1))-1;n++){ MT3[k][l][m]+=MT1[k][n][l]*MT2[m][n];} } } }
              }

           if(i==3 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[k][l][0]*MT2[0][m];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=MT1.Dimension(2)-1;n++){ MT3[k][l][m]+=MT1[k][l][n]*MT2[n][m];} } } }
              }
           if(i==3 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ MT3[k][l][m]=MT1[k][l][0]*MT2[m][0];} } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=1;n<=MT1.Dimension(2)-1;n++){ MT3[k][l][m]+=MT1[k][l][n]*MT2[m][n];} } } }
              }

           return MT3;
          }

//partial specialization for second rank tensor and fourth rank
template <class Type>
TARRAY<Type,4> Contract(const TARRAY<Type,2> &MT1,const TARRAY<Type,4> &MT2,int i,int j)
         { if(MT1.Dimension(i-1)!=MT2.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,2>&,const TARRAY<Type,4>&,int,int)");}

           std::vector<size_t> dimensions3(MT1.Dimensions());
           dimensions3.erase(dimensions3.begin()+i-1);
           for(int k=0;k<=static_cast<int>(MT2.Rank())-1;k++){ dimensions3.push_back(MT2.Dimension(k));}
           dimensions3.erase(dimensions3.begin()+j);    // j=(2-1) + (j-1)

           TARRAY<Type,4> MT3(dimensions3);

           if(i==1 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[0][k]*MT2[0][l][m][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(0))-1;p++){ MT3[k][l][m][n]+=MT1[p][k]*MT2[p][l][m][n];} } } } }
              }
           if(i==1 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[0][k]*MT2[l][0][m][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(0))-1;p++){ MT3[k][l][m][n]+=MT1[p][k]*MT2[l][p][m][n];} } } } }
              }
           if(i==1 && j==3)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[0][k]*MT2[l][m][0][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(0))-1;p++){ MT3[k][l][m][n]+=MT1[p][k]*MT2[l][m][p][n];} } } } }
              }
           if(i==1 && j==4)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[0][k]*MT2[l][m][n][0];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(0))-1;p++){ MT3[k][l][m][n]+=MT1[p][k]*MT2[l][m][n][p];} } } } }
              }

           if(i==2 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][0]*MT2[0][l][m][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(1))-1;p++){ MT3[k][l][m][n]+=MT1[k][p]*MT2[p][l][m][n];} } } } }
              }
           if(i==2 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][0]*MT2[l][0][m][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(1))-1;p++){ MT3[k][l][m][n]+=MT1[k][p]*MT2[l][p][m][n];} } } } }
              }
           if(i==2 && j==3)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][0]*MT2[l][m][0][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(1))-1;p++){ MT3[k][l][m][n]+=MT1[k][p]*MT2[l][m][p][n];} } } } }
              }
           if(i==2 && j==4)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][0]*MT2[l][m][n][0];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(1))-1;p++){ MT3[k][l][m][n]+=MT1[k][p]*MT2[l][m][n][p];} } } } }
              }

           return MT3;
          }

//partial specialization for fourth rank tensor and second rank
template <class Type>
TARRAY<Type,4> Contract(const TARRAY<Type,4> &MT1,const TARRAY<Type,2> &MT2,int i,int j)
         { if(MT1.Dimension(i-1)!=MT2.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,2>&,const TARRAY<Type,4>&,int,int)");}

           std::vector<size_t> dimensions3(MT1.Dimensions());
           dimensions3.erase(dimensions3.begin()+i-1);

           for(int k=0;k<=static_cast<int>(MT2.Rank())-1;k++){ dimensions3.push_back(MT2.Dimension(k));}
           dimensions3.erase(dimensions3.begin()+2+j);  // 2+j=(4-1) + (j-1)

           TARRAY<Type,4> MT3(dimensions3);

           if(i==1 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[0][k][l][m]*MT2[0][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(0))-1;p++){ MT3[k][l][m][n]+=MT1[p][k][l][m]*MT2[p][n];} } } } }
              }
           if(i==2 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][0][l][m]*MT2[0][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(1))-1;p++){ MT3[k][l][m][n]+=MT1[k][p][l][m]*MT2[p][n];} } } } }
              }
           if(i==3 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][l][0][m]*MT2[0][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=MT1.Dimension(2)-1;p++){ MT3[k][l][m][n]+=MT1[k][l][p][m]*MT2[p][n];} } } } }
              }
           if(i==4 && j==1)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][l][m][0]*MT2[0][n];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=MT1.Dimension(3)-1;p++){ MT3[k][l][m][n]+=MT1[k][l][m][p]*MT2[p][n];} } } } }
              }

           if(i==1 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[0][k][l][m]*MT2[n][0];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(0))-1;p++){ MT3[k][l][m][n]+=MT1[p][k][l][m]*MT2[n][p];} } } } }
              }
           if(i==2 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][0][l][m]*MT2[n][0];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=static_cast<int>(MT1.Dimension(1))-1;p++){ MT3[k][l][m][n]+=MT1[k][p][l][m]*MT2[n][p];} } } } }
              }
           if(i==3 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][l][0][m]*MT2[n][0];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=MT1.Dimension(2)-1;p++){ MT3[k][l][m][n]+=MT1[k][l][p][m]*MT2[n][p];} } } } }
              }
           if(i==4 && j==2)
             { for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ MT3[k][l][m][n]=MT1[k][l][m][0]*MT2[n][0];} } } }
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=0;l<=static_cast<int>(MT3.Dimension(1))-1;l++){ for(int m=0;m<=static_cast<int>(MT3.Dimension(2))-1;m++){ for(int n=0;n<=MT3.Dimension(3)-1;n++){ for(int p=1;p<=MT1.Dimension(3)-1;p++){ MT3[k][l][m][n]+=MT1[k][l][m][p]*MT2[n][p];} } } } }
              }

           return MT3;
          }

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

// These versions, all involve first rank tensors, don't acccept the first rank tensor's index
template <class Type,size_t rank1>
TARRAY<Type,rank1+1-2> Contract(const TARRAY<Type,rank1> &MT1,const TARRAY<Type,1> &MT2,int i)
         { if(MT1.Dimension(i-1)!=MT2.Dimension()){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type>&,const TARRAY<Type,1>&,int)");}

           if(i==1)
             { std::vector<size_t> dimensions3(MT1.Dimensions()); dimensions3.erase(dimensions3.begin());
               TARRAY<Type,rank1+1-2> MT3(dimensions3);

               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ MT3[k]=MT1[0][k]*MT2[0];}
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=1;l<=static_cast<int>(MT1.Dimension(0))-1;l++){ MT3[k]+=MT1[l][k]*MT2[l];} }
               return MT3;
              }
           else{ TARRAY<Type,rank1> MT3(MT1);
                 MT3.SwapIndices(i,1);
                 return Contract(MT3,MT2,1);
                }
          }

//partial specialization
template <class Type>
TARRAY<Type,1> Contract(const TARRAY<Type,2> &MT1,const TARRAY<Type,1> &MT2,int i)
         { if(MT1.Dimension(i-1)!=MT2.Dimension()){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,2>&,const TARRAY<Type,1>&,int)");}

           if(i==1)
             { TARRAY<Type,1> MT3(MT2.Dimension(1));
               for(int k=0;k<=static_cast<int>(MT3.Dimension())-1;k++){ MT3[k]=MT1[0][k]*MT2[0];}
               for(int k=0;k<=static_cast<int>(MT3.Dimension())-1;k++){ for(int l=1;l<=static_cast<int>(MT1.Dimension(0))-1;l++){ MT3[k]+=MT1[l][k]*MT2[l];} }
               return MT3;
              }
           else{ TARRAY<Type,2> MT3(MT1);
                 MT3.SwapIndices(i,1);
                 return Contract(MT3,MT2,1);
                }
          }

//**********************************

template <class Type,size_t rank2>
TARRAY<Type,1+rank2-2> Contract(const TARRAY<Type,1> &MT1,const TARRAY<Type,rank2> &MT2,int j)
         { if(MT1.Dimension()!=MT2.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,1>&,const TARRAY<Type>&,int)");}

           if(j==1)
             { std::vector<size_t> dimensions3(MT2.Dimensions()); dimensions3.erase(dimensions3.begin());
               TARRAY<Type,1+rank2-2> MT3(dimensions3);
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ MT3[k]=MT1[0]*MT2[0][k];}
               for(int k=0;k<=static_cast<int>(MT3.Dimension(0))-1;k++){ for(int l=1;l<=static_cast<int>(MT1.Dimension(0))-1;l++){ MT3[k]+=MT1[l]*MT2[l][k];} }
               return MT3;
              }
           else{ TARRAY<Type,rank2> MT3(MT2);
                 MT3.SwapIndices(j,1);
                 return Contract(MT1,MT3,1);
                }
          }

//partial specialization
template <class Type>
TARRAY<Type,1> Contract(const TARRAY<Type,1> &MT1,const TARRAY<Type,2> &MT2,int j)
         { if(MT1.Dimension()!=MT2.Dimension(j-1)){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,1>&,const TARRAY<Type,2>&,int)");}

           if(j==1)
             { TARRAY<Type,1> MT3(MT2.Dimension(1));
               for(int k=0;k<=static_cast<int>(MT3.Dimension())-1;k++){ MT3[k]=MT1[0]*MT2[0][k];}
               for(int k=0;k<=static_cast<int>(MT3.Dimension())-1;k++){ for(int l=1;l<=static_cast<int>(MT1.Dimension())-1;l++){ MT3[k]+=MT1[l]*MT2[l][k];} }
               return MT3;
              }
           else{ TARRAY<Type,1> MT3(MT2.Dimension(0));
                 for(int k=0;k<=static_cast<int>(MT3.Dimension())-1;k++){ MT3[k]=MT1[0]*MT2[k][0];}
                 for(int k=0;k<=static_cast<int>(MT3.Dimension())-1;k++){ for(int l=1;l<=static_cast<int>(MT1.Dimension())-1;l++){ MT3[k]+=MT1[l]*MT2[k][l];} }
                 return MT3;
                }
          }

//**********************************

template <class Type>
TARRAY<Type,0> Contract(const TARRAY<Type,1> &MT1,const TARRAY<Type,1> &MT2)
         { if(MT1.Dimension()!=MT2.Dimension()){ throw DIFFERENT_DIMENSIONS("Contract(const TARRAY<Type,1>&,const TARRAY<Type,1>&)");}
           TARRAY<Type,0> MT3(MT1[0]*MT2[0]);
           for(int k=1;k<=static_cast<int>(MT1.Dimension())-1;k++){ MT3+=MT1[k]*MT2[k];}
           return MT3;
          }

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

template <class Type>
TARRAY<Type,2> Invert(TARRAY<Type,2> MT)
         { MATRIX<Type> M(MT.Dimension(0),MT.Dimension(1));
           for(int i=0;i<=static_cast<int>(MT.Dimension(0))-1;i++){ for(int j=0;j<=static_cast<int>(MT.Dimension(1))-1;j++){ M[i][j]=MT[i][j];} }
           M.LUInvert();
           for(int i=0;i<=static_cast<int>(MT.Dimension(0))-1;i++){ for(int j=0;j<=static_cast<int>(MT.Dimension(1))-1;j++){ MT[i][j]=M[i][j];} }
           return MT;
          }

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

template <class Type>
std::ostream& operator<<(std::ostream &os,const TARRAY<Type,3> &MT)
         { os<<"\n";
           for(int i=0;i<=static_cast<int>(MT.Dimension(0))-2;i++)
              { os<<MT[i]<<"\n----";}
           os<<MT[MT.Dimension(0)-1];
           return os;
          }

template <class Type>
std::ostream& operator<<(std::ostream &os,const TARRAY<Type,2> &MT)
         { for(int i=0;i<=static_cast<int>(MT.Dimension(0))-1;i++)
              { os<<"\n";
                for(int j=0;j<=static_cast<int>(MT.Dimension(1))-1;j++){ os<<MT[i][j]<<"\t";}
               }
           return os;
          }

template <class Type>
std::ostream& operator<<(std::ostream &os,const TARRAY<Type,1> &MT)
         { os<<"\n"; for(int i=0;i<=static_cast<int>(MT.Dimension())-1;i++){ os<<MT[i]<<"\t";} return os;}

template <class Type>
std::ostream& operator<<(std::ostream &os,const TARRAY<Type,0> &MT){ os<<(Type)MT; return os;}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

#endif
