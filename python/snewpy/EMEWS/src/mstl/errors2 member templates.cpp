#if !defined(_BASIC_ERRORS_MEMBERS)
#define _BASIC_ERRORS_MEMBERS

template <class Type> void EQUAL_VALUES<Type>::Message(void)
     { std::cout<<std::endl<<BASIC_ERROR::what()<<" "<<which()<<" "<<who()<<" "<<where();
       std::cout<<"\n\tNumber "<<Number()<<std::flush;
      }

// **************************************************************************

template <class Type> void NOT_EQUAL_VALUES<Type>::Message(void)
     { std::cout<<std::endl<<BASIC_ERROR::what()<<" "<<which()<<" "<<who()<<" "<<where();
       std::cout<<"\n\tFirst Number "<<FirstNumber()<<"\tSecond Number "<<SecondNumber()<<std::flush;
      }

// **************************************************************************

template <class Type> void POLE<Type>::Message(void)
     { std::cout<<std::endl<<BASIC_ERROR::what()<<" "<<which()<<" "<<who()<<" "<<where();
       std::cout<<"\n\tNumber "<<Number()<<std::flush;
      }

// **************************************************************************

template <class Type1,class Type2> void INCOMPATABLE_TYPES<Type1,Type2>::Message(void)
     { std::cout<<std::endl<<BASIC_ERROR::what()<<" "<<which()<<" "<<who()<<" "<<where();
       std::cout<<" "<<typeid(Type1).name()<<" "<<typeid(Type2).name()<<std::flush;
      }

// **************************************************************************

template <class Type> void DID_NOT_CONVERGE<Type>::Message(void)
     { std::cout<<std::endl<<BASIC_ERROR::what()<<" "<<which()<<" "<<who()<<" "<<where();
       std::cout<<"\n\tValue Reached is "/*<<ValueReached()*/<<" after "<<Iterations()<<" iterations."<<std::flush;
      }

// **************************************************************************

template <class Type1,class Type2> void CONVERGED_TO<Type1,Type2>::Message(void)
     { std::cout<<std::endl<<BASIC_ERROR::what()<<" "<<which()<<" "<<who()<<" "<<where();
       std::cout<<"\n\tValue Reached is "/*<<ValueReached()*/<<" after "<<Interval()<<" interval."<<std::flush;
      }

// **************************************************************************

template <class Type> void OUT_OF_RANGE<Type>::Message(void)
     { std::cout<<std::endl<<BASIC_ERROR::what()<<" "<<which()<<" "<<who()<<" "<<where();
       std::cout<<"\n\t: The "<<typeid(Type).name()<<" "<<Value()<<" is not between "<<LowerLimit()<<" and "<<UpperLimit()<<std::flush;
      }

#endif


