
#include <functional>
#include <vector>

#include "mstl.h"

#if !defined(_FUNCTIONAL2)
#define _FUNCTIONAL2

// this a structure similar to unary_function and binary_function in the std and is 
// designed to be inherited for nullary functions (functions that take no arguements)
template <typename RType> struct nullary_function;
template <typename AType1,typename AType2,typename AType3,typename RType> struct ternary_function;
template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType> struct quaternary_function;

// this a structure that is similar to unary_function and binary_function
// in the std but is meant for classes and member functions
template <class CType,typename RType> struct nullary_member_function;
template <class CType,typename AType,typename RType> struct unary_member_function;
template <class CType,typename AType1,typename AType2,typename RType> struct binary_member_function;
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> struct ternary_member_function;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> struct quaternary_member_function;

// *******************************************************************

// similar to the std classes 
template <typename RType> class pointer_to_nullary_function;
template <typename AType1,typename AType2,typename AType3,typename RType> class pointer_to_ternary_function;
template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class pointer_to_quaternary_function;

// overloading of ptr_fun
template <typename RType> inline pointer_to_nullary_function<RType> ptr_fun(RType(*fptr)(void));
template <typename AType1,typename AType2,typename AType3,typename RType> inline pointer_to_ternary_function<AType1,AType2,AType3,RType> ptr_fun(RType(*fptr)(AType1,AType2,AType3));
template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType> inline pointer_to_quaternary_function<AType1,AType2,AType3,AType4,RType> ptr_fun(RType(*fptr)(AType1,AType2,AType3,AType4));

// similar to the std classes but for member functions
template <class CType,typename RType> class pointer_to_nullary_member_function;
template <class CType,typename AType,typename RType> class pointer_to_unary_member_function;
template <class CType,typename AType1,typename AType2,typename RType> class pointer_to_binary_member_function;
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> class pointer_to_ternary_member_function;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class pointer_to_quaternary_member_function;

// *******************************************************************
// *******************************************************************
// *******************************************************************

// additional unary function classes
template <typename Type> struct equal;
template <typename Type> struct plus;
template <typename Type> struct square;
template <typename Type> struct cube;

// additional binary function classes
template <typename Type> struct exclusive_or;

// *******************************************************************
// *******************************************************************
// ****** FUNCTION AND MEMBER FUNCTION BASE CLASSES ******************
// *******************************************************************
// *******************************************************************

// adds a virtual operator() to nullary_function and the std unary_function and binary_function in functional
template <typename RType> class NULLARYFUNCTORBASE;
template <typename AType,typename RType=AType> class UNARYFUNCTORBASE;
template <typename AType1,typename AType2=AType1,typename RType=AType1> class BINARYFUNCTORBASE;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1> class TERNARYFUNCTORBASE;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> class QUATERNARYFUNCTORBASE;

// *******************************************************************
// *******************************************************************
// ******** FUNCTION AND MEMBER FUNCTION CLASSES *********************
// *******************************************************************
// *******************************************************************

// turns a functor or a function into a functor
template <typename RType,typename FType> class NULLARYFUNCTOR;
template <typename AType,typename RType,typename FType> class UNARYFUNCTOR;
template <typename AType1,typename AType2,typename RType,typename FType> class BINARYFUNCTOR;
template <typename AType1,typename AType2,typename AType3,typename RType,typename FType> class TERNARYFUNCTOR;
template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> class QUATERNARYFUNCTOR;

// *******************************************************************
// *******************************************************************
// **** BOUND FUNCTION AND MEMBER FUNCTION CLASSES *******************
// *******************************************************************
// *******************************************************************

// turns a unary functor into a nullary functor by binding the arguement
template <typename AType,typename RType=AType,typename FType=RType(*)(AType)> class UNARY2NULLARYFUNCTOR;

// turns a binary functor into a nullary functor by binding both arguements
template <typename AType1,typename AType2=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2)> class BINARY2NULLARYFUNCTOR;
// turns a binary functor into a nullary functor by binding one of the arguements
template <typename AType1,typename AType2=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2)> class BINARY2UNARYFUNCTOR1st;
template <typename AType1,typename AType2=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2)> class BINARY2UNARYFUNCTOR2nd;

// turns a ternary functor into a nullary functor by binding all arguements
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3)> class TERNARY2NULLARYFUNCTOR;
// turns a ternary functor into a unary functor by binding two of the arguements
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3)> class TERNARY2UNARYFUNCTOR1st2nd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3)> class TERNARY2UNARYFUNCTOR2nd3rd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3)> class TERNARY2UNARYFUNCTOR1st3rd;
// turns a ternary functor into a binary functor by binding one of the arguements
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3)> class TERNARY2BINARYFUNCTOR1st;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3)> class TERNARY2BINARYFUNCTOR2nd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3)> class TERNARY2BINARYFUNCTOR3rd;

// turns a quaternary functor into a nullary functor by binding all arguements
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2NULLARYFUNCTOR;
// turns a quaternary functor into a unary functor by binding three of the arguements
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2UNARYFUNCTOR1st2nd3rd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2UNARYFUNCTOR1st2nd4th;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2UNARYFUNCTOR1st3rd4th;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2UNARYFUNCTOR2nd3rd4th;
// turns a quaternary functor into a binary functor by binding two of the arguements
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2BINARYFUNCTOR1st2nd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2BINARYFUNCTOR1st3rd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2BINARYFUNCTOR1st4th;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2BINARYFUNCTOR2nd3rd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2BINARYFUNCTOR2nd4th;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2BINARYFUNCTOR3rd4th;
// turns a quaternary functor into a binary functor by binding one of the arguements
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2TERNARYFUNCTOR1st;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2TERNARYFUNCTOR2nd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2TERNARYFUNCTOR3rd;
template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)> class QUATERNARY2TERNARYFUNCTOR4th;

// turns a unary member functor into a nullary functor by binding the arguement
template <class CType,typename AType,typename RType> class UNARY2NULLARYMEMBERFUNCTION;
// turns a binary member functor into a nullary functor by binding both arguements
template <class CType,typename AType1,typename AType2,typename RType> class BINARY2NULLARYMEMBERFUNCTION;
// turns a binary member functor into a nullary functor by binding one of the arguements
template <class CType,typename AType1,typename AType2,typename RType> class BINARY2UNARYMEMBERFUNCTION1st;
template <class CType,typename AType1,typename AType2,typename RType> class BINARY2UNARYMEMBERFUNCTION2nd;
// turns a ternary functor into a nullary functor by binding all arguements
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> class TERNARY2NULLARYMEMBERFUNCTION;
// turns a ternary functor into a unary functor by binding two of the arguements
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> class TERNARY2UNARYMEMBERFUNCTION1st2nd;
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> class TERNARY2UNARYMEMBERFUNCTION2nd3rd;
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> class TERNARY2UNARYMEMBERFUNCTION1st3rd;
// turns a ternary functor into a binary functor by binding one of the arguements
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> class TERNARY2BINARYMEMBERFUNCTION1st;
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> class TERNARY2BINARYMEMBERFUNCTION2nd;
template <class CType,typename AType1,typename AType2,typename AType3,typename RType> class TERNARY2BINARYMEMBERFUNCTION3rd;

// turns a quaternary functor into a nullary functor by binding all arguements
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2NULLARYMEMBERFUNCTION;
// turns a quaternary functor into a unary functor by binding three of the arguements
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2UNARYMEMBERFUNCTION1st2nd3rd;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2UNARYMEMBERFUNCTION1st2nd4th;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2UNARYMEMBERFUNCTION1st3rd4th;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2UNARYMEMBERFUNCTION2nd3rd4th;
// turns a quaternary functor into a binary functor by binding two of the arguements
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2BINARYMEMBERFUNCTION1st2nd;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2BINARYMEMBERFUNCTION1st3rd;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2BINARYMEMBERFUNCTION1st4th;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2BINARYMEMBERFUNCTION2nd3rd;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2BINARYMEMBERFUNCTION2nd4th;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2BINARYMEMBERFUNCTION3rd4th;
// turns a quaternary functor into a binary functor by binding one of the arguements
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2TERNARYMEMBERFUNCTION1st;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2TERNARYMEMBERFUNCTION2nd;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2TERNARYMEMBERFUNCTION3rd;
template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class QUATERNARY2TERNARYMEMBERFUNCTION4th;

// *******************************************************************
// *******************************************************************
// *******************************************************************

// this a structure similar to unary_function and binary_function in the std and is 
// designed to be inherited for nullary functions (functions that take no arguements)
template <typename RType> struct nullary_function { typedef RType result_type;};

template <typename AType1,typename AType2,typename AType3,typename RType> struct ternary_function 
         { typedef AType1 arguement_type1;
           typedef AType2 arguement_type2;
           typedef AType3 arguement_type3;
           typedef RType result_type;
          };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType> struct quaternary_function 
         { typedef AType1 arguement_type1;
           typedef AType2 arguement_type2;
           typedef AType3 arguement_type3;
           typedef AType4 arguement_type4;
           typedef RType result_type;
          };

// *******************************************************************
// *******************************************************************
// *******************************************************************

// this a structure that is similar to unary_function and binary_function
// in the std but is meant for classes and member functions
template <class CType,typename RType> struct nullary_member_function : public nullary_function<RType> 
          { typedef CType class_type;};

template <class CType,typename AType,typename RType> struct unary_member_function : public std::unary_function<AType,RType>
         { typedef CType class_type;};

template <class CType,typename AType1,typename AType2,typename RType> struct binary_member_function : public std::binary_function<AType1,AType2,RType>
         { typedef CType class_type;};

template <class CType,typename AType1,typename AType2,typename AType3,typename RType> struct ternary_member_function : public ternary_function<AType1,AType2,AType3,RType>
         { typedef CType class_type;};

template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> struct quaternary_member_function : public quaternary_function<AType1,AType2,AType3,AType4,RType>
         { typedef CType class_type;};

// *******************************************************************
// *******************************************************************
// *******************************************************************

// overloading of ptr_fun
template <typename RType> inline pointer_to_nullary_function<RType> ptr_fun(RType(*fptr)(void))
         { return pointer_to_nullary_function<RType>(fptr);}

template <typename AType1,typename AType2,typename AType3,typename RType> inline pointer_to_ternary_function<AType1,AType2,AType3,RType> ptr_fun(RType(*fptr)(AType1,AType2,AType3))
         { return pointer_to_ternary_function<AType1,AType2,AType3,RType>(fptr);}

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType> inline pointer_to_quaternary_function<AType1,AType2,AType3,AType4,RType> ptr_fun(RType(*fptr)(AType1,AType2,AType3,AType4))
         { return pointer_to_quaternary_function<AType1,AType2,AType3,AType4,RType>(fptr);}

// *******************************************************************
// *******************************************************************
// *******************************************************************

template <typename RType> class pointer_to_nullary_function : public nullary_function<RType>
         { protected : RType(*f)(void);

           public: pointer_to_nullary_function(void){;}
                   explicit pointer_to_nullary_function(RType(*F)(void)) : f(F) {;}

                   RType operator()(void) const { return (*f)(); }
          };

template <typename AType1,typename AType2,typename AType3,typename RType> class pointer_to_ternary_function : public ternary_function<AType1,AType2,AType3,RType>
         { protected : RType(*f)(AType1,AType2,AType3);

           public: pointer_to_ternary_function(void){;}
                   explicit pointer_to_ternary_function(RType(*F)(AType1,AType2,AType3)) : f(F) {;}

                   RType operator()(AType1 A1,AType2 A2,AType3 A3) const { return (*f)(A1,A2,A3); }
          };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType> class pointer_to_quaternary_function : public quaternary_function<AType1,AType2,AType3,AType4,RType>
         { protected : RType(*f)(AType1,AType2,AType3,AType4);

           public: pointer_to_quaternary_function(void){;}
                   explicit pointer_to_quaternary_function(RType(*F)(AType1,AType2,AType3,AType4)) : f(F) {;}

                   RType operator()(AType1 A1,AType2 A2,AType3 A3,AType4 A4) const { return (*f)(A1,A2,A3,A4); }
          };

// *******************************************************
// *******************************************************
// *******************************************************

template <class CType,typename RType> 
class pointer_to_nullary_member_function : public nullary_member_function<CType,RType>
         { protected : RType(CType::*mfptr)(void);
                       const CType &c;

           public: pointer_to_nullary_member_function(void){ mfptr=NULL;}   
                   explicit pointer_to_nullary_member_function(const CType &C,RType(CType::*MFPTR)(void)) 
                            : c(C), mfptr(MFPTR) {;}

                   RType operator()(void) const { return (c.*mfptr)();}
          };

// *******************************************************

template <class CType,typename AType,typename RType>
class pointer_to_unary_member_function : public unary_member_function<CType,AType,RType>
         { private : RType(CType::*mfptr)(AType);
                     const CType &c;
                     
           public  : pointer_to_unary_member_function(void){ mfptr=NULL;}     
                     explicit pointer_to_unary_member_function(const CType &C,RType(CType::*MFPTR)(AType)) 
                              : c(C), mfptr(MFPTR) {;}

                     RType operator()(AType ARGS) const { return (c.*mfptr)(ARGS);}
          };

// *******************************************************

template <class CType,typename AType1,typename AType2,typename RType>
class pointer_to_binary_member_function : public binary_member_function<CType,AType1,AType2,RType>
         { private : RType(CType::*mfptr)(AType1,AType2);
                     const CType &c;
                     
           public  : pointer_to_binary_member_function(void){ mfptr=NULL;}
                     explicit pointer_to_binary_member_function(const CType &C,RType(CType::*MFPTR)(AType1,AType2))
                              : c(C), mfptr(MFPTR) {;}                     

                     RType operator()(AType1 A1,AType2 A2) const { return (c.*mfptr)(A1,A2);}
          };

// *******************************************************

template <class CType,typename AType1,typename AType2,typename AType3,typename RType> 
class pointer_to_ternary_member_function : public ternary_member_function<CType,AType1,AType2,AType3,RType>
         { protected : RType(CType::*mfptr)(AType1,AType2,AType3);
                       const CType &c;

           public: pointer_to_ternary_member_function(void){ mfptr=NULL;}   
                   explicit pointer_to_ternary_member_function(const CType &C,RType(CType::*MFPTR)(AType1,AType2,AType3)) 
                            : c(C), mfptr(MFPTR) {;}

                   RType operator()(AType1 A1,AType2 A2,AType3 A3) const { return (c.*mfptr)(A1,A2,A3);}
          };

// *******************************************************

template <class CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType> 
class pointer_to_quaternary_member_function : public quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
         { protected : RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
                       const CType &c;

           public: pointer_to_quaternary_member_function(void){ mfptr=NULL;}   
                   explicit pointer_to_quaternary_member_function(const CType &C,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                            : c(C), mfptr(MFPTR) {;}

                   RType operator()(AType1 A1,AType2 A2,AType3 A3,AType4 A4) const { return (c.*mfptr)(A1,A2,A3,A4);}
          };

// *******************************************************
// *******************************************************
// *******************************************************

template <typename Type> struct equal : public std::unary_function<Type,Type>{
         Type operator()(const Type &x) const { return x;}
        };

template <typename Type> struct plus : public std::unary_function<Type,Type>{
         Type operator()(const Type &x) const { return x;}
        };

template <typename Type> struct square : public std::unary_function<Type,Type>{
         Type operator()(const Type &x) const { return x*x;}
        };

template <typename Type> struct cube : public std::unary_function<Type,Type>{
         Type operator()(const Type &x) const { return x*x*x;}
        };

// *******************************************************
// *******************************************************
// *******************************************************

template <typename Type> struct exclusive_or : public std::binary_function<Type,Type,Type>{
         Type operator()(const Type &x,const Type &y) const { return x^y;}
        };

// *******************************************************
// *******************************************************
// *******************************************************

template <typename RType>
class NULLARYFUNCTORBASE : 
public nullary_function<RType>
         { public  : virtual ~NULLARYFUNCTORBASE(void){;}
                     virtual RType operator()(void) const =0;
          };

template <typename AType,typename RType>
class UNARYFUNCTORBASE : 
public std::unary_function<AType,RType>
         { public  : virtual ~UNARYFUNCTORBASE(void){;}
                     virtual RType operator()(AType) const =0;
          };

template <typename AType1,typename AType2,typename RType>
class BINARYFUNCTORBASE : 
public std::binary_function<AType1,AType2,RType>
         { public  : virtual ~BINARYFUNCTORBASE(void){;}
                     virtual RType operator()(AType1,AType2) const =0;
          };

template <typename AType1,typename AType2,typename AType3,typename RType>
class TERNARYFUNCTORBASE : 
public ternary_function<AType1,AType2,AType3,RType>
         { public  : virtual ~TERNARYFUNCTORBASE(void){;}
                     virtual RType operator()(AType1,AType2,AType3) const =0;
          };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType>
class QUATERNARYFUNCTORBASE : 
public quaternary_function<AType1,AType2,AType3,AType4,RType>
         { public  : virtual ~QUATERNARYFUNCTORBASE(void){;}
                     virtual RType operator()(AType1,AType2,AType3,AType4) const =0;
          };

// *******************************************************
// *******************************************************
// *******************************************************

template <typename CType,typename RType>
class NULLARYMEMBERFUNCTIONBASE : public nullary_member_function<CType,RType>
         { public  : virtual ~NULLARYMEMBERFUNCTIONBASE(void){;}
                     virtual RType operator()(void) const =0;
          };

template <typename CType,typename AType,typename RType>
class UNARYMEMBERFUNCTIONBASE : public unary_member_function<CType,AType,RType>
         { public  : virtual ~UNARYMEMBERFUNCTIONBASE(void){;}
                     virtual RType operator()(AType) const =0;
          };

template <typename CType,typename AType1,typename AType2,typename RType>
class BINARYMEMBERFUNCTIONBASE : public binary_member_function<CType,AType1,AType2,RType>
         { public  : virtual ~BINARYMEMBERFUNCTIONBASE(void){;}
                     virtual RType operator()(AType1,AType2) const =0;
          };

template <typename CType,typename AType1,typename AType2,typename AType3,typename RType>
class TERNARYMEMBERFUNCTIONBASE : public ternary_member_function<CType,AType1,AType2,AType3,RType>
         { public  : virtual ~TERNARYMEMBERFUNCTIONBASE(void){;}
                     virtual RType operator()(AType1,AType2,AType3) const =0;
          };

template <typename CType,typename AType1,typename AType2,typename AType3,typename AType4,typename RType>
class QUATERNARYMEMBERFUNCTIONBASE : public quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
         { public  : virtual ~QUATERNARYMEMBERFUNCTIONBASE(void){;}
                     virtual RType operator()(AType1,AType2,AType3,AType4) const =0;
          };

// *******************************************************
// *******************************************************
// *******************************************************
// *******************************************************
// *******************************************************
// *******************************************************

template <typename RType,typename FType=RType(*)(void)>
class NULLARYFUNCTOR : public virtual NULLARYFUNCTORBASE<RType>
         { private : const FType &f;
                     
           public  : NULLARYFUNCTOR(const FType &F) : f(F) {;}
                     virtual ~NULLARYFUNCTOR(void){;}

                     virtual RType operator()(void) const { return f();}
          };

// *******************************************************
// *******************************************************
// *******************************************************

template <typename AType,typename RType=AType,typename FType=RType(*)(AType)>
class UNARYFUNCTOR : public virtual UNARYFUNCTORBASE<AType,RType>
         { private : const FType &f;
                     
           public  : UNARYFUNCTOR(const FType &F) : f(F) {;}
                     virtual ~UNARYFUNCTOR(void){;}

                     virtual RType operator()(AType A) const { return f(A);}
          };

// *******************************************************
// *******************************************************
// *******************************************************

template <typename AType1,typename AType2=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2)>
class BINARYFUNCTOR : public virtual BINARYFUNCTORBASE<AType1,AType2,RType>
         { private : const FType &f;
                     
           public  : BINARYFUNCTOR(const FType &F) : f(F) {;}
                     virtual ~BINARYFUNCTOR(void){;}

                     virtual RType operator()(AType1 A1,AType2 A2) const { return f(A1,A2);}
          };

// *******************************************************
// *******************************************************
// *******************************************************

template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3)>
class TERNARYFUNCTOR : public virtual TERNARYFUNCTORBASE<AType1,AType2,AType3,RType>
         { private : const FType &f;
                     
           public  : TERNARYFUNCTOR(const FType &F) : f(F) {;}
                     virtual ~TERNARYFUNCTOR(void){;}

                     virtual RType operator()(AType1 A1,AType2 A2,AType3 A3) const { return f(A1,A2,A3);}
          };

// *******************************************************
// *******************************************************
// *******************************************************

template <typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1,typename FType=RType(*)(AType1,AType2,AType3,AType4)>
class QUATERNARYFUNCTOR : public virtual QUATERNARYFUNCTORBASE<AType1,AType2,AType3,AType4,RType>
         { private : const FType &f;
                     
           public  : QUATERNARYFUNCTOR(const FType &F) : f(F) {;}
                     virtual ~QUATERNARYFUNCTOR(void){;}

                     virtual RType operator()(AType1 A1,AType2 A2,AType3 A3,AType4 A4) const { return f(A1,A2,A3,A4);}
          };
         
// *******************************************************
// *******************************************************
// *******************************************************

// turns a unary functor into a nullary functor by binding the arguement
template <typename AType,typename RType,typename FType> 
class UNARY2NULLARYFUNCTOR : public virtual NULLARYFUNCTORBASE<RType>
      { private : const AType &a;                  
                  const FType &f;

        public : UNARY2NULLARYFUNCTOR(const FType &F,const AType &A) : a(A), f(F) {;}
                 ~UNARY2NULLARYFUNCTOR(void){;}
                 
                 void ChangeBoundArguement(const AType &A){ a=A;}
                 RType operator()(void) const { return f(a);}
       };
          
// turns a binary functor into a nullary functor by binding both arguements
template <typename AType1,typename AType2,typename RType,typename FType> 
class BINARY2NULLARYFUNCTOR : public virtual NULLARYFUNCTORBASE<RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const FType &f;
      
        public : BINARY2NULLARYFUNCTOR(const FType &F,const AType1 &A1,const AType2 &A2) : a1(A1), a2(A2), f(F) {;}
                 ~BINARY2NULLARYFUNCTOR(void){;}

                 void ChangeBoundArguement1st(const AType1 &A1){ a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2){ a2=A2;}
                 
                 RType operator()(void) const { return f(a1,a2);}
       };

// turns a binary functor into a unary functor by binding the first arguement       
template <typename AType1,typename AType2,typename RType,typename FType> 
class BINARY2UNARYFUNCTOR1st : public virtual UNARYFUNCTORBASE<AType2,RType>
      { private : const AType1 &a1;
                  const FType &f;
      
        public : BINARY2UNARYFUNCTOR1st(const FType &F,const AType1 &A1) : a1(A1), f(F) {;}
                 ~BINARY2UNARYFUNCTOR1st(void){;}

                 void ChangeBoundArguement1st(const AType1 &A1){ a1=A1;}
                 
                 RType operator()(AType2 A2) const { return f(a1,A2);}
       };

// turns a binary functor into a unary functor by binding the second arguement              
template <typename AType1,typename AType2,typename RType,typename FType> 
class BINARY2UNARYFUNCTOR2nd : public virtual UNARYFUNCTORBASE<AType1,RType>
      { private : const AType2 &a2;
                  const FType &f;
      
        public : BINARY2UNARYFUNCTOR2nd(const FType &F,const AType2 &A2) : a2(A2), f(F) {;}
                 ~BINARY2UNARYFUNCTOR2nd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType1 A1) const { return f(A1,a2);}
       };

// turns a ternary functor into a nullary functor by binding all arguements
template <typename AType1,typename AType2,typename AType3,typename RType,typename FType> 
class TERNARY2NULLARYFUNCTOR : public virtual NULLARYFUNCTORBASE<RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const AType3 &a3;
                  const FType &f;
      
        public : TERNARY2NULLARYFUNCTOR(const FType &F,const AType1 &A1,const AType2 &A2,const AType3 &A3) : a1(A1), a2(A2), a3(A3), f(F) {;}
                 ~TERNARY2NULLARYFUNCTOR(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(void) const { return f(a1,a2,a3);}
       };

// turns a ternary functor into a unary functor by binding two of the arguements
template <typename AType1,typename AType2,typename AType3,typename RType,typename FType> 
class TERNARY2UNARYFUNCTOR1st2nd : public virtual UNARYFUNCTORBASE<AType3,RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const FType &f;
      
        public : TERNARY2UNARYFUNCTOR1st2nd(const FType &F,const AType1 &A1,const AType2 &A2) : a1(A1), a2(A2), f(F) {;}
                 ~TERNARY2UNARYFUNCTOR1st2nd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType3 A3) const { return f(a1,a2,A3);}
       };

template <typename AType1,typename AType2,typename AType3,typename RType,typename FType> 
class TERNARY2UNARYFUNCTOR2nd3rd : public virtual UNARYFUNCTORBASE<AType1,RType>
      { private : const AType2 &a2; 
                  const AType3 &a3;
                  const FType &f;
      
        public : TERNARY2UNARYFUNCTOR2nd3rd(const FType &F,const AType2 &A2,const AType3 &A3) : a2(A2), a3(A3), f(F) {;}
                 ~TERNARY2UNARYFUNCTOR2nd3rd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType1 A1) const { return f(A1,a2,a3);}
       };

template <typename AType1,typename AType2,typename AType3,typename RType,typename FType> 
class TERNARY2UNARYFUNCTOR1st3rd : public virtual UNARYFUNCTORBASE<AType2,RType>
      { private : const AType1 &a1; 
                  const AType3 &a3;
                  const FType &f;
      
        public : TERNARY2UNARYFUNCTOR1st3rd(const FType &F,const AType1 &A1,const AType3 &A3) : a1(A1), a3(A3), f(F) {;}
                 ~TERNARY2UNARYFUNCTOR1st3rd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType2 A2) const { return f(a1,A2,a3);}
       };

template <typename AType1,typename AType2,typename AType3,typename RType,typename FType> 
class TERNARY2BINARYFUNCTOR1st : public virtual BINARYFUNCTORBASE<AType2,AType3,RType>
      { private : const AType1 &a1; 
                  const FType &f;
      
        public : TERNARY2BINARYFUNCTOR1st(const FType &F,const AType1 &A1) : a1(A1), f(F) {;}
                 ~TERNARY2BINARYFUNCTOR1st(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 
                 RType operator()(AType2 A2,AType3 A3) const { return f(a1,A2,A3);}
       };

template <typename AType1,typename AType2,typename AType3,typename RType,typename FType> 
class TERNARY2BINARYFUNCTOR2nd : public virtual BINARYFUNCTORBASE<AType1,AType3,RType>
      { private : const AType2 &a2; 
                  const FType &f;
      
        public : TERNARY2BINARYFUNCTOR2nd(const FType &F,const AType2 &A2) : a2(A2), f(F) {;}
                 ~TERNARY2BINARYFUNCTOR2nd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType1 A1,AType3 A3) const { return f(A1,a2,A3);}
       };

template <typename AType1,typename AType2,typename AType3,typename RType,typename FType> 
class TERNARY2BINARYFUNCTOR3rd : public virtual BINARYFUNCTORBASE<AType1,AType2,RType>
      { private : const AType1 &a3;
                  const FType &f; 
      
        public : TERNARY2BINARYFUNCTOR3rd(const FType &F,const AType3 &A3) : a3(A3), f(F) {;}
                 ~TERNARY2BINARYFUNCTOR3rd(void) {;}

                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType1 A1,AType2 A2) const { return f(A1,A2,a3);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2NULLARYFUNCTOR : public virtual NULLARYFUNCTORBASE<RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const AType3 &a3;
                  const AType4 &a4;
                  const FType &f;
      
        public : QUATERNARY2NULLARYFUNCTOR(const FType &F,const AType1 &A1,const AType2 &A2,const AType3 &A3,const AType4 &A4) : a1(A1), a2(A2), a3(A3), a4(A4), f(F) {;}
                 ~QUATERNARY2NULLARYFUNCTOR(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(void) const { return f(a1,a2,a3,a4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2UNARYFUNCTOR1st2nd3rd : public virtual UNARYFUNCTORBASE<AType4,RType> 
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const AType3 &a3;
                  const FType &f;
      
        public : QUATERNARY2UNARYFUNCTOR1st2nd3rd(const FType &F,const AType1 &A1,const AType2 &A2,const AType3 &A3) : a1(A1), a2(A2), a3(A3), f(F) {;}
                 ~QUATERNARY2UNARYFUNCTOR1st2nd3rd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType4 A4) const { return f(a1,a2,a3,A4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2UNARYFUNCTOR1st2nd4th : public virtual UNARYFUNCTORBASE<AType3,RType> 
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const AType4 &a4;
                  const FType &f;
      
        public : QUATERNARY2UNARYFUNCTOR1st2nd4th(const FType &F,const AType1 &A1,const AType2 &A2,const AType4 &A4) : a1(A1), a2(A2), a4(A4), f(F) {;}
                 ~QUATERNARY2UNARYFUNCTOR1st2nd4th(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType3 A3) const { return f(a1,a2,A3,a4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2UNARYFUNCTOR1st3rd4th : public virtual UNARYFUNCTORBASE<AType2,RType> 
      { private : const AType1 &a1; 
                  const AType3 &a3;
                  const AType4 &a4;
                  const FType &f;
      
        public : QUATERNARY2UNARYFUNCTOR1st3rd4th(const FType &F,const AType1 &A1,const AType3 &A3,const AType4 &A4) : a1(A1), a3(A3), a4(A4), f(F) {;}
                 ~QUATERNARY2UNARYFUNCTOR1st3rd4th(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType2 A2) const { return f(a1,A2,a3,a4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2UNARYFUNCTOR2nd3rd4th : public virtual UNARYFUNCTORBASE<AType1,RType> 
      { private : const AType2 &a2;
                  const AType3 &a3;
                  const AType4 &a4;
                  const FType &f;
      
        public : QUATERNARY2UNARYFUNCTOR2nd3rd4th(const FType &F,const AType2 &A2,const AType3 &A3,const AType4 &A4) : a2(A2), a3(A3), a4(A4), f(F) {;}
                 ~QUATERNARY2UNARYFUNCTOR2nd3rd4th(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType1 A1) const { return f(A1,a2,a3,a4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2BINARYFUNCTOR1st2nd : public virtual BINARYFUNCTORBASE<AType3,AType4,RType> 
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const FType &f;
      
        public : QUATERNARY2BINARYFUNCTOR1st2nd(const FType &F,const AType1 &A1,const AType2 &A2) : a1(A1), a2(A2), f(F) {;}
                 ~QUATERNARY2BINARYFUNCTOR1st2nd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType3 A3,AType4 A4) const { return f(a1,a2,A3,A4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2BINARYFUNCTOR1st3rd : public virtual BINARYFUNCTORBASE<AType2,AType4,RType> 
      { private : const AType1 &a1; 
                  const AType3 &a3;
                  const FType &f;
      
        public : QUATERNARY2BINARYFUNCTOR1st3rd(const FType &F,const AType1 &A1,const AType3 &A3) : a1(A1), a3(A3), f(F) {;}
                 ~QUATERNARY2BINARYFUNCTOR1st3rd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType2 A2,AType4 A4) const { return f(a1,A2,a3,A4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2BINARYFUNCTOR1st4th : public virtual BINARYFUNCTORBASE<AType2,AType3,RType> 
      { private : const AType1 &a1; 
                  const AType4 &a4;
                  const FType &f;
      
        public : QUATERNARY2BINARYFUNCTOR1st4th(const FType &F,const AType1 &A1,const AType4 &A4) : a1(A1), a4(A4), f(F) {;}
                 ~QUATERNARY2BINARYFUNCTOR1st4th(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType2 A2,AType3 A3) const { return f(a1,A2,A3,a4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2BINARYFUNCTOR2nd3rd : public virtual BINARYFUNCTORBASE<AType1,AType4,RType> 
      { private : const AType2 &a2;
                  const AType3 &a3;
                  const FType &f;
      
        public : QUATERNARY2BINARYFUNCTOR2nd3rd(const FType &F,const AType2 &A2,const AType3 &A3) : a2(A2), a3(A3), f(F) {;}
                 ~QUATERNARY2BINARYFUNCTOR2nd3rd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType1 A1,AType4 A4) const { return f(A1,a2,a3,A4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2BINARYFUNCTOR2nd4th : public virtual BINARYFUNCTORBASE<AType1,AType3,RType> 
      { private : const AType2 &a2;
                  const AType4 &a4;
                  const FType &f;
      
        public : QUATERNARY2BINARYFUNCTOR2nd4th(const FType &F,const AType2 &A2,const AType4 &A4) : a2(A2), a4(A4), f(F) {;}
                 ~QUATERNARY2BINARYFUNCTOR2nd4th(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType1 A1,AType3 A3) const { return f(A1,a2,A3,a4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2BINARYFUNCTOR3rd4th : public virtual BINARYFUNCTORBASE<AType1,AType2,RType> 
      { private : const AType3 &a3;
                  const AType4 &a4;
                  const FType &f;
      
        public : QUATERNARY2BINARYFUNCTOR3rd4th(const FType &F,const AType3 &A3,const AType4 &A4) : a3(A3), a4(A4), f(F) {;}
                 ~QUATERNARY2BINARYFUNCTOR3rd4th(void) {;}

                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType1 A1,AType2 A2) const { return f(A1,A2,a3,a4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2TERNARYFUNCTOR1st : public virtual TERNARYFUNCTORBASE<AType2,AType3,AType4,RType> 
      { private : const AType1 &a1; 
                  const FType &f;
      
        public : QUATERNARY2TERNARYFUNCTOR1st(const FType &F,const AType1 &A1) : a1(A1), f(F) {;}
                 ~QUATERNARY2TERNARYFUNCTOR1st(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 
                 RType operator()(AType2 A2,AType3 A3,AType4 A4) const { return f(a1,A2,A3,A4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2TERNARYFUNCTOR2nd : public virtual TERNARYFUNCTORBASE<AType1,AType3,AType4,RType>
      { private : const AType2 &a2;
                  const FType &f;
      
        public : QUATERNARY2TERNARYFUNCTOR2nd(const FType &F,const AType2 &A2) : a2(A2), f(F) {;}
                 ~QUATERNARY2TERNARYFUNCTOR2nd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType1 A1,AType3 A3,AType4 A4) const { return f(A1,a2,A3,A4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2TERNARYFUNCTOR3rd : public virtual TERNARYFUNCTORBASE<AType1,AType2,AType4,RType>
      { private : const AType3 &a3;
                  const FType &f;
      
        public : QUATERNARY2TERNARYFUNCTOR3rd(const FType &F,const AType3 &A3) : a3(A3), f(F) {;}
                 ~QUATERNARY2TERNARYFUNCTOR3rd(void) {;}

                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType1 A1,AType2 A2,AType4 A4) const { return f(A1,A2,a3,A4);}
       };

template <typename AType1,typename AType2,typename AType3,typename AType4,typename RType,typename FType> 
class QUATERNARY2TERNARYFUNCTOR4th : public virtual TERNARYFUNCTORBASE<AType1,AType2,AType3,RType>
      { private : const AType4 &a4;
                  const FType &f;
      
        public : QUATERNARY2TERNARYFUNCTOR4th(const FType &F,const AType4 &A4) : a4(A4), f(F) {;}
                 ~QUATERNARY2TERNARYFUNCTOR4th(void) {;}

                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType1 A1,AType2 A2,AType3 A3) const { return f(A1,A2,A3,a4);}
       };
       
// *******************************************************
// *******************************************************
// *******************************************************

// turns a unary member functor into a nullary functor by binding the arguement
template <class CType,typename AType,typename RType=AType> 
class UNARY2NULLARYMEMBERFUNCTOR : public virtual unary_member_function<CType,AType,RType>
      { private : const AType &a;
                  const CType &c;
                  RType(CType::*mfptr)(AType);
      
        public : UNARY2NULLARYMEMBERFUNCTOR(const CType &C,const AType &A,RType(CType::*MFPTR)(AType)) 
                   : a(A), c(C), mfptr(MFPTR) {;}//{ mfptr=MFPTR;}
                 ~UNARY2NULLARYMEMBERFUNCTOR(void) {;}

                 void ChangeBoundArguement(const AType &A) { a=A;}
                 
                 RType operator()(void) const { return (c.*mfptr)(a);}
       };   
          
// turns a binary member functor into a nullary functor by binding both arguements
template <class CType,typename AType1,typename AType2=AType1,typename RType=AType1> 
class BINARY2NULLARYMEMBERFUNCTOR : public virtual binary_member_function<CType,AType1,AType2,RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2);
      
        public : BINARY2NULLARYMEMBERFUNCTOR(const CType &C,const AType1 &A1,const AType2 &A2,RType(CType::*MFPTR)(AType1,AType2)) 
                    : a1(A1), a2(A2), c(C) { mfptr=MFPTR;}
                 ~BINARY2NULLARYMEMBERFUNCTOR(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(void) const { return (c.*mfptr)(a1,a2);}
       };
       
// turns a binary member functor into a unary functor by binding the first arguement       
template <class CType,typename AType1,typename AType2=AType1,typename RType=AType1> 
class BINARY2UNARYMEMBERFUNCTOR1st : public virtual binary_member_function<CType,AType1,AType2,RType>
      { private : const AType1 &a1;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2);
      
        public : BINARY2UNARYMEMBERFUNCTOR1st(const CType &C,const AType1 &A1,RType(CType::*MFPTR)(AType1,AType2)) 
                    : a1(A1), c(C) { mfptr=MFPTR;}
                 ~BINARY2UNARYMEMBERFUNCTOR1st(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 
                 RType operator()(AType2 A2) const { return (c.*mfptr)(a1,A2);}
       };

// turns a binary member functor into a unary functor by binding the second arguement              
template <class CType,typename AType1,typename AType2=AType1,typename RType=AType1> 
class BINARY2UNARYMEMBERFUNCTOR2nd : public virtual binary_member_function<CType,AType1,AType2,RType>
      { private : const AType2 &a2;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2);
      
        public : BINARY2UNARYMEMBERFUNCTOR2nd(const CType &C,const AType2 &A2,RType(CType::*MFPTR)(AType1,AType2))
                    : a2(A2), c(C) { mfptr=MFPTR;}
                 ~BINARY2UNARYMEMBERFUNCTOR2nd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType1 A1) const { return (c.*mfptr)(A1,a2);}
       };      

// turns a ternary functor into a nullary functor by binding all arguements
template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1> 
class TERNARY2NULLARYMEMBERFUNCTOR : public virtual ternary_member_function<CType,AType1,AType2,AType3,RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const AType3 &a3;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3);
      
        public : TERNARY2NULLARYMEMBERFUNCTOR(const CType &C,const AType1 &A1,const AType2 &A2,const AType3 &A3,RType(CType::*MFPTR)(AType1,AType2,AType3)) 
                    : a1(A1), a2(A2), a3(A3), c(C) { mfptr=MFPTR;}
                 ~TERNARY2NULLARYMEMBERFUNCTOR(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(void) const { return (c.*mfptr)(a1,a2,a3);}
       };

// turns a ternary functor into a unary functor by binding two of the arguements
template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1> 
class TERNARY2UNARYMEMBERFUNCTOR1st2nd : public virtual ternary_member_function<CType,AType1,AType2,AType3,RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3);
      
        public : TERNARY2UNARYMEMBERFUNCTOR1st2nd(const CType &C,const AType1 &A1,const AType2 &A2,RType(CType::*MFPTR)(AType1,AType2,AType3)) 
                    : a1(A1), a2(A2), c(C) { mfptr=MFPTR;}
                 ~TERNARY2UNARYMEMBERFUNCTOR1st2nd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType3 A3) const { return (c.*mfptr)(a1,a2,A3);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1> 
class TERNARY2UNARYMEMBERFUNCTOR2nd3rd : public virtual ternary_member_function<CType,AType1,AType2,AType3,RType>
      { private : const AType2 &a2; 
                  const AType3 &a3;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3);
      
        public : TERNARY2UNARYMEMBERFUNCTOR2nd3rd(const CType &C,const AType2 &A2,const AType3 &A3,RType(CType::*MFPTR)(AType1,AType2,AType3)) 
                    : a2(A2), a3(A3), c(C) { mfptr=MFPTR;}
                 ~TERNARY2UNARYMEMBERFUNCTOR2nd3rd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType1 A1) const { return (c.*mfptr)(A1,a2,a3);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1> 
class TERNARY2UNARYMEMBERFUNCTOR1st3rd : public virtual ternary_member_function<CType,AType1,AType2,AType3,RType>
      { private : const AType1 &a1; 
                  const AType3 &a3;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3);
      
        public : TERNARY2UNARYMEMBERFUNCTOR1st3rd(const CType &C,const AType1 &A1,const AType3 &A3,RType(CType::*MFPTR)(AType1,AType2,AType3)) 
                    : a1(A1), a3(A3), c(C) { mfptr=MFPTR;}
                 ~TERNARY2UNARYMEMBERFUNCTOR1st3rd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType2 A2) const { return (c.*mfptr)(a1,A2,a3);}
       };

// turns a ternary functor into a binary functor by binding one of the arguements
template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1> 
class TERNARY2BINARYMEMBERFUNCTOR1st : public virtual ternary_member_function<CType,AType1,AType2,AType3,RType>
      { private : const AType1 &a1; 
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3);
      
        public : TERNARY2BINARYMEMBERFUNCTOR1st(const CType &C,const AType1 &A1,RType(CType::*MFPTR)(AType1,AType2,AType3)) 
                    : a1(A1), c(C) { mfptr=MFPTR;}
                 ~TERNARY2BINARYMEMBERFUNCTOR1st(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 
                 RType operator()(AType2 A2,AType3 A3) const { return (c.*mfptr)(a1,A2,A3);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1> 
class TERNARY2BINARYMEMBERFUNCTOR2nd : public virtual ternary_member_function<CType,AType1,AType2,AType3,RType>
      { private : const AType2 &a2; 
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3);
      
        public : TERNARY2BINARYMEMBERFUNCTOR2nd(const CType &C,const AType2 &A2,RType(CType::*MFPTR)(AType1,AType2,AType3)) 
                    : a2(A2), c(C) { mfptr=MFPTR;}
                 ~TERNARY2BINARYMEMBERFUNCTOR2nd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType1 A1,AType3 A3) const { return (c.*mfptr)(A1,a2,A3);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename RType=AType1> 
class TERNARY2BINARYMEMBERFUNCTOR3rd : public virtual ternary_member_function<CType,AType1,AType2,AType3,RType>
      { private : const AType3 &a3; 
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3);
      
        public : TERNARY2BINARYMEMBERFUNCTOR3rd(const CType &C,const AType3 &A3,RType(CType::*MFPTR)(AType1,AType2,AType3)) 
                    : a3(A3), c(C) { mfptr=MFPTR;}
                 ~TERNARY2BINARYMEMBERFUNCTOR3rd(void) {;}

                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType1 A1,AType2 A2) const { return (c.*mfptr)(A1,A2,a3);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2NULLARYMEMBERFUNCTOR : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const AType3 &a3;
                  const AType4 &a4;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2NULLARYMEMBERFUNCTOR(const CType &C,const AType1 &A1,const AType2 &A2,const AType3 &A3,const AType4 &A4,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a1(A1), a2(A2), a3(A3), a4(A4), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2NULLARYMEMBERFUNCTOR(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(void) const { return (c.*mfptr)(a1,a2,a3,a4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2UNARYMEMBERFUNCTOR1st2nd3rd : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const AType3 &a3;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2UNARYMEMBERFUNCTOR1st2nd3rd(const CType &C,const AType1 &A1,const AType2 &A2,const AType3 &A3,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a1(A1), a2(A2), a3(A3), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2UNARYMEMBERFUNCTOR1st2nd3rd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType4 A4) const { return (c.*mfptr)(a1,a2,a3,A4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2UNARYMEMBERFUNCTOR1st2nd4th : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const AType4 &a4;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2UNARYMEMBERFUNCTOR1st2nd4th(const CType &C,const AType1 &A1,const AType2 &A2,const AType4 &A4,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a1(A1), a2(A2), a4(A4) { mfptr=MFPTR;}
                 ~QUATERNARY2UNARYMEMBERFUNCTOR1st2nd4th(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType3 A3) const { return (c.*mfptr)(a1,a2,A3,a4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2UNARYMEMBERFUNCTOR1st3rd4th : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType1 &a1; 
                  const AType3 &a3;
                  const AType4 &a4;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2UNARYMEMBERFUNCTOR1st3rd4th(const CType &C,const AType1 &A1,const AType3 &A3,const AType4 &A4,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a1(A1), a3(A3), a4(A4), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2UNARYMEMBERFUNCTOR1st3rd4th(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType2 A2) const { return (c.*mfptr)(a1,A2,a3,a4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2UNARYMEMBERFUNCTOR2nd3rd4th : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType2 &a2;
                  const AType3 &a3;
                  const AType4 &a4;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2UNARYMEMBERFUNCTOR2nd3rd4th(const CType &C,const AType2 &A2,const AType3 &A3,const AType4 &A4,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a2(A2), a3(A3), a4(A4), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2UNARYMEMBERFUNCTOR2nd3rd4th(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType1 A1) const { return (c.*mfptr)(A1,a2,a3,a4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2BINARYMEMBERFUNCTOR1st2nd : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType1 &a1; 
                  const AType2 &a2;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2BINARYMEMBERFUNCTOR1st2nd(const CType &C,const AType1 &A1,const AType2 &A2,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a1(A1), a2(A2), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2BINARYMEMBERFUNCTOR1st2nd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType3 A3,AType4 A4) const { return (c.*mfptr)(a1,a2,A3,A4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2BINARYMEMBERFUNCTOR1st3rd : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType1 &a1; 
                  const AType3 &a3;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2BINARYMEMBERFUNCTOR1st3rd(const CType &C,const AType1 &A1,const AType3 &A3,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a1(A1), a3(A3), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2BINARYMEMBERFUNCTOR1st3rd(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType2 A2,AType4 A4) const { return (c.*mfptr)(a1,A2,a3,A4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2BINARYMEMBERFUNCTOR1st4th : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType1 &a1; 
                  const AType4 &a4;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2BINARYMEMBERFUNCTOR1st4th(const CType &C,const AType1 &A1,const AType4 &A4,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a1(A1), a4(A4), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2BINARYMEMBERFUNCTOR1st4th(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType2 A2,AType3 A3) const { return (c.*mfptr)(a1,A2,A3,a4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2BINARYMEMBERFUNCTOR2nd3rd : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType2 &a2;
                  const AType3 &a3;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2BINARYMEMBERFUNCTOR2nd3rd(const CType &C,const AType2 &A2,const AType3 &A3,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a2(A2), a3(A3), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2BINARYMEMBERFUNCTOR2nd3rd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType1 A1,AType4 A4) const { return (c.*mfptr)(A1,a2,a3,A4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2BINARYMEMBERFUNCTOR2nd4th : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType2 &a2;
                  const AType4 &a4;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2BINARYMEMBERFUNCTOR2nd4th(const CType &C,const AType2 &A2,const AType4 &A4,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a2(A2), a4(A4), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2BINARYMEMBERFUNCTOR2nd4th(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType1 A1,AType3 A3) const { return (c.*mfptr)(A1,a2,A3,a4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2BINARYMEMBERFUNCTOR3rd4th : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType3 &a3;
                  const AType4 &a4;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2BINARYMEMBERFUNCTOR3rd4th(const CType &C,const AType3 &A3,const AType4 &A4,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a3(A3), a4(A4), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2BINARYMEMBERFUNCTOR3rd4th(void) {;}

                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType1 A1,AType2 A2) const { return (c.*mfptr)(A1,A2,a3,a4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2TERNARYMEMBERFUNCTOR1st : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType1 &a1; 
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2TERNARYMEMBERFUNCTOR1st(const CType &C,const AType1 &A1,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a1(A1), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2TERNARYMEMBERFUNCTOR1st(void) {;}

                 void ChangeBoundArguement1st(const AType1 &A1) { a1=A1;}
                 
                 RType operator()(AType2 A2,AType3 A3,AType4 A4) const { return (c.*mfptr)(a1,A2,A3,A4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2TERNARYMEMBERFUNCTOR2nd : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType2 &a2;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2TERNARYMEMBERFUNCTOR2nd(const CType &C,const AType2 &A2,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a2(A2), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2TERNARYMEMBERFUNCTOR2nd(void) {;}

                 void ChangeBoundArguement2nd(const AType2 &A2) { a2=A2;}
                 
                 RType operator()(AType1 A1,AType3 A3,AType4 A4) const { return (c.*mfptr)(A1,a2,A3,A4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2TERNARYMEMBERFUNCTOR3rd : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType3 &a3;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);
      
        public : QUATERNARY2TERNARYMEMBERFUNCTOR3rd(const CType &C,const AType3 &A3,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a3(A3), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2TERNARYMEMBERFUNCTOR3rd(void) {;}

                 void ChangeBoundArguement3rd(const AType3 &A3) { a3=A3;}
                 
                 RType operator()(AType1 A1,AType2 A2,AType4 A4) const { return (c.*mfptr)(A1,A2,a3,A4);}
       };

template <class CType,typename AType1,typename AType2=AType1,typename AType3=AType1,typename AType4=AType1,typename RType=AType1> 
class QUATERNARY2TERNARYMEMBERFUNCTOR4th : public virtual quaternary_member_function<CType,AType1,AType2,AType3,AType4,RType>
      { private : const AType4 &a4;
                  const CType &c;
                  RType(CType::*mfptr)(AType1,AType2,AType3,AType4);

      
        public : QUATERNARY2TERNARYMEMBERFUNCTOR4th(const CType &C,const AType4 &A4,RType(CType::*MFPTR)(AType1,AType2,AType3,AType4)) 
                    : a4(A4), c(C) { mfptr=MFPTR;}
                 ~QUATERNARY2TERNARYMEMBERFUNCTOR4th(void) {;}

                 void ChangeBoundArguement4th(const AType4 &A4) { a4=A4;}
                 
                 RType operator()(AType1 A1,AType2 A2,AType3 A3) const { return (c.*mfptr)(A1,A2,A3,a4);}
       };

// *******************************************************
// *******************************************************
// *******************************************************

#include "mstl.h"

#endif

