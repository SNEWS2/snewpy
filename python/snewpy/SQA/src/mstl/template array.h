
#include<cmath>
#include<cstdlib>
#include<cstdarg>
#include<iostream>
#include<typeinfo>
#include<vector>

#include "mstl.h"

#if !defined(_TARRAY)
#define _TARRAY

// *******************************************************************
// ************************* TEMPLATE ARRAY **************************
// *******************************************************************

template <class Type,int rank> class TARRAY;

// ******************

// Outer product
template <class Type,int rank1,int rank2>
TARRAY<Type,rank1+rank2> operator*(const TARRAY<Type,rank1>&,const TARRAY<Type,rank2>&);

template <class Type,int rank>
TARRAY<Type,rank> operator*(const TARRAY<Type,rank>&,const Type&);
template <class Type,int rank>
TARRAY<Type,rank> operator*(const Type&,const TARRAY<Type,rank>&);

// ******************

//Addition and Subtraction
template <class Type,int rank>
TARRAY<Type,rank> operator+(const TARRAY<Type,rank>&,const TARRAY<Type,rank>&);

template <class Type,int rank>
TARRAY<Type,rank> operator-(const TARRAY<Type,rank>&,const TARRAY<Type,rank>&);

// ******************

template <class Type,int rank>
TARRAY<Type,rank> SwapAdjacentIndices(const TARRAY<Type,rank>&,int);

// *******************************

// Contraction of indicii on one tensor
template <class Type,int rank> TARRAY<Type,rank-2> Contract(const TARRAY<Type,rank>&,int,int);

// Contraction of indicii on two tensors
template <class Type,int rank1,int rank2> TARRAY<Type,rank1+rank2-2> Contract(const TARRAY<Type,rank1>&,const TARRAY<Type,rank2>&,int,int);
//partial specializations
template <class Type> TARRAY<Type,0> Contract(const TARRAY<Type,1> &MT1,const TARRAY<Type,1> &MT2,int,int);
template <class Type,int rank2> TARRAY<Type,1+rank2-2> Contract(const TARRAY<Type,1> &MT1,const TARRAY<Type,rank2> &MT2,int,int j);
template <class Type,int rank1> TARRAY<Type,rank1+1-2> Contract(const TARRAY<Type,rank1> &MT1,const TARRAY<Type,1> &MT2,int i,int);
template <class Type> TARRAY<Type,2> Contract(const TARRAY<Type,2> &MT1,const TARRAY<Type,2> &MT2,int i,int j);
template <class Type> TARRAY<Type,3> Contract(const TARRAY<Type,2> &MT1,const TARRAY<Type,3> &MT2,int i,int j);
template <class Type> TARRAY<Type,3> Contract(const TARRAY<Type,3> &MT1,const TARRAY<Type,2> &MT2,int i,int j);
template <class Type> TARRAY<Type,4> Contract(const TARRAY<Type,2> &MT1,const TARRAY<Type,4> &MT2,int i,int j);
template <class Type> TARRAY<Type,4> Contract(const TARRAY<Type,4> &MT1,const TARRAY<Type,2> &MT2,int i,int j);

// These three allow contraction of a tensor with a first rank tensor without specifying
// the index for the first rank (its 1 after all)
template <class Type,int rank1> TARRAY<Type,rank1+1-2> Contract(const TARRAY<Type,rank1>&,const TARRAY<Type,1>&,int);
template <class Type,int rank2> TARRAY<Type,1+rank2-2> Contract(const TARRAY<Type,1>&,const TARRAY<Type,rank2>&,int);
template <class Type> TARRAY<Type,0> Contract(const TARRAY<Type,1>&,const TARRAY<Type,1>&);

// *************** Functions for specific tensors ********************

template <class Type> TARRAY<Type,2> Invert(TARRAY<Type,2> MT);

template <class Type> std::ostream& operator<<(std::ostream&,const TARRAY<Type,3>&);
template <class Type> std::ostream& operator<<(std::ostream&,const TARRAY<Type,2>&);
template <class Type> std::ostream& operator<<(std::ostream&,const TARRAY<Type,1>&);
template <class Type> std::ostream& operator<<(std::ostream&,const TARRAY<Type,0>&);

// *******************************************************************
// *******************************************************************
// ************************* MTENSOR *********************************
// *******************************************************************
// *******************************************************************

template <class Type,int rank> class TARRAY
      	{ private :
           bool empty;
           std::vector<int> dimensions;
           TARRAY<Type,rank-1> *mt;

           void Create(void){ empty=true; mt=NULL; dimensions=std::vector<int>(rank);}
           template <int rank2> void Create(TARRAY<Type,rank2>*,std::vector<int> DIMENSIONS); // create from a list of dimensions
           void Create(TARRAY<Type,2>*,std::vector<int> DIMENSIONS);

           void Copy(Type *T);
           void Copy(const TARRAY<Type,rank>&);

           void Fill(const Type &T);

           template <int rank2> void Destroy(TARRAY<Type,rank2>*);
           void Destroy(TARRAY<Type,2>*);

           // Swap index i with index i+1
           template <int rank2> TARRAY<Type,rank>& SwapAdjacentIndices(TARRAY<Type,rank2>*,int);
           TARRAY<Type,2>& SwapAdjacentIndices(TARRAY<Type,2>*,int);

           TARRAY<Type,rank-1>& Reference(int i) { return mt[i];}
           TARRAY<Type,rank-1> Value(int i) const { return mt[i];}

           template <int rank2> Type& Reference(TARRAY<Type,rank2>*,const std::vector<int> &i,int i0) { return Reference(&mt[i[i0]],i,i0+1);}
           template <int rank2> Type Value(TARRAY<Type,rank2>*,const std::vector<int> &i,int i0) const { return Value(&mt[i[i0]],i,i0+1);}

           public :
           TARRAY(void){ Create();}
           TARRAY(const TARRAY<Type,rank> &MT);
           explicit TARRAY(std::vector<int> DIMENSIONS);
           explicit TARRAY(std::vector<int> DIMENSIONS,const Type &T);

           ~TARRAY(void){ Destroy(this);}

           bool Empty(void) const { return empty;}
           int Rank(void) const { return rank;}
           int Dimension(int i) const;
           std::vector<int> Dimensions(void) const { return dimensions;}
           std::vector<int> size(void) const { return dimensions;}

           // iterators
           typedef Type*       iterator;
           typedef const Type* const_iterator;

           iterator begin(void){ return mt;}             const_iterator begin(void) const { return mt;}
           iterator end(void){ return mt+Dimension();}   const_iterator end(void) const { return mt+Dimension();}

           TARRAY<Type,rank>& SwapIndices(int,int);     // Swap index i with index j
           TARRAY<Type,rank>& SwapAdjacentIndices(int); // Swap index i with index i+1

           TARRAY<Type,rank>& operator=(const TARRAY<Type,rank> &MT);

           TARRAY<Type,rank>& operator+=(const TARRAY<Type,rank>&);
           TARRAY<Type,rank>& operator-=(const TARRAY<Type,rank>&);

           TARRAY<Type,rank>& operator*=(const Type&);
           TARRAY<Type,rank>& operator/=(const Type&);

           TARRAY<Type,rank-1>& operator[](const int&);
           TARRAY<Type,rank-1> operator[](const int&) const;

           Type& operator[](const std::vector<int>&);
           Type operator[](const std::vector<int>&) const;

           friend class TARRAY<Type,rank+1>;
       	 };

// *******************************************************************
// *******************************************************************
// *******************************************************************

template <class Type> class TARRAY<Type,1>
      	{ private : bool empty;
                    int dimension;
                    Type *mt;

           void Create(void){ empty=true; dimension=0; mt=NULL;}
           void Create(const int &);

           void Copy(Type*);
           void Copy(const std::vector<Type>&);
           void Copy(const TARRAY<Type,1>&);

           void Fill(const Type &T);

           void Destroy(void);

           Type& Reference(int i) { return mt[i];}
           Type Value(int i) const { return mt[i];}

           Type& Reference(TARRAY<Type,1>*,std::vector<int> &i,int i0) { return mt[i[i0]];}
           Type Value(TARRAY<Type,1>*,std::vector<int> &i,int i0) const { return mt[i[i0]];}

           public:
           TARRAY(void){ Create();}
           explicit TARRAY(int DIMENSION){ Create(DIMENSION);}
           explicit TARRAY(int DIMENSION,const Type &T){ Create(DIMENSION); Fill(T);}
           explicit TARRAY(int DIMENSION,Type *T){ Create(DIMENSION); Copy(T);}
           explicit TARRAY(std::vector<int> DIMENSIONS);
           TARRAY(const TARRAY<Type,1> &MT);

           ~TARRAY(void){ Destroy();}

           bool Empty(void) const { return empty;}
           int Rank(void) const { return 1;}
           int Dimension(void) const;
           std::vector<int> Dimensions(void) const { return std::vector<int>(1,dimension);}

           TARRAY<Type,1>& operator=(const TARRAY<Type,1> &MT);
           TARRAY<Type,1>& operator=(const std::vector<Type> &v);

           TARRAY<Type,1>& operator+=(const TARRAY<Type,1>&);
           TARRAY<Type,1>& operator-=(const TARRAY<Type,1>&);

           TARRAY<Type,1>& operator*=(const Type&);
           TARRAY<Type,1>& operator/=(const Type&);

           Type& operator[](const int&);
           Type operator[](const int&) const;

           Type& operator[](const std::vector<int>&);
           Type operator[](const std::vector<int>&) const;

           friend class TARRAY<Type,2>;
          };

// *******************************************************************
// *******************************************************************
// *******************************************************************

template <class Type> class TARRAY<Type,0>
      	{ private : Type t;

           void Copy(const TARRAY<Type,0> &MT){ t=MT.t;}

           public:
           TARRAY(void){;}
           TARRAY(Type T){ t=T;}
           TARRAY(const TARRAY<Type,0> &MT){ Copy(MT);}

           int Rank(void) const { return 0;}
           int Dimension(void) const { return 1;}

           operator Type() const { return t;}

           TARRAY<Type,0>& operator=(const TARRAY<Type,0> &MT){ t=MT.t; return (*this);}
           TARRAY<Type,0>& operator=(Type T){ t=T; return (*this);}

           TARRAY<Type,0>& operator+=(const TARRAY<Type,0> &MT){ t+=MT.t; return (*this);}
           TARRAY<Type,0>& operator-=(const TARRAY<Type,0> &MT){ t-=MT.t; return (*this);}

           TARRAY<Type,0>& operator+=(Type T){ t+=T; return (*this);}
           TARRAY<Type,0>& operator-=(Type T){ t-=T; return (*this);}
           TARRAY<Type,0>& operator*=(Type T){ t*=T; return (*this);}
           TARRAY<Type,0>& operator/=(Type T){ t/=T; return (*this);}
          };

#endif

