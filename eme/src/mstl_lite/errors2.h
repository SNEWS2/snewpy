
#include<cstdio>
#include<stdexcept>
#include<typeinfo>
#include<iostream>
#include<string>
#include<iomanip>

#include "mstl.h"

#if !defined(_BASIC_ERRORS)
#define _BASIC_ERRORS

// ***************************************************************************************************
// ***************************************************************************************************
// ***************************************************************************************************

class BASIC_ERROR;

// ****************************************************************************

class CANNOT_FIND;

class CANNOT_READ;

class CANNOT_WRITE;

// ****************************************************************************

class EMPTY;

class DIFFERENT_DIMENSIONS;

class DIFFERENT_LENGTHS;

class DIFFERENT_SIZES;

// ****************************************************************************

class DIVISION_BY_ZERO;

class EVEN_NUMBER;

class NEGATIVE_INTEGER;

class NEGATIVE_NUMBER;

class NOT_AN_INTEGER;

class ODD_NUMBER;

class ZERO_NUMBER;

class FUNCTION_ERROR;

class NO_SOLUTION;

// ****************************************************************************

class UNKNOWN_ERROR;

// ****************************************************************************

template <class Type> class EQUAL_VALUES;

template <class Type> class NOT_EQUAL_VALUES;

template <class Type> class POLE;

template <class Type1,class Type2> class INCOMPATABLE_TYPES;

template <class Type> class DID_NOT_CONVERGE;

template <class Type1,class Type2> class CONVERGED_TO;

// **************************************************************************

template <class Type> class OUT_OF_RANGE;

// **************************************************************************
// **************************************************************************
// **************************************************************************

class BASIC_ERROR : public std::exception
      { private :
        std::string function, object, line;

        protected :
        static std::string Unknown,ThrownIn,InObject,AtLine;

        BASIC_ERROR(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
        : std::exception(){ function=ThrownIn+WHICH;
                            object=InObject+WHO;
                            line=AtLine+WHERE;
                           }
        BASIC_ERROR(const std::string &WHICH,const std::string &WHO,const int &WHERE)
        : std::exception(){ function=ThrownIn+WHICH;
                            object=InObject+WHO;
                            char buffer[7]; sprintf(buffer,"%d",WHERE);
                            line=AtLine+std::string(buffer);
                           }

        public :
        BASIC_ERROR(const BASIC_ERROR &BE)
        : std::exception(BE){ function=BE.function; object=BE.object; line=BE.line;}
        
        virtual ~BASIC_ERROR(void) throw() {;}

        void ChangeWhich(const std::string &WHICH){ function=ThrownIn+WHICH;}
        void ChangeFunction(const std::string &WHICH){ function=ThrownIn+WHICH;}

        void ChangeWho(const std::string &WHO){ object=InObject+WHO;}
        void ChangeObject(const std::string &WHO){ object=InObject+WHO;}

        void ChangeWhere(const std::string &WHERE){ line=AtLine+WHERE;}
        void ChangeLine(const std::string &WHERE){ line=AtLine+WHERE;}
        void ChangeWhere(const int &WHERE){ char buffer[7]; sprintf(buffer,"%d",WHERE); line=AtLine+std::string(buffer);}
        void ChangeLine(const int &WHERE){ char buffer[7]; sprintf(buffer,"%d",WHERE); line=AtLine+std::string(buffer);}

        void Change(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
                 { function=( WHICH==Unknown ? function : ThrownIn+WHICH );
                   object=( WHO==Unknown ? object : InObject+WHO );
                   line=( WHERE==Unknown ? line : AtLine+WHERE );
                  }
        void Change(const std::string &WHICH,const std::string &WHO,const int &WHERE)
                 { function=ThrownIn+WHICH;
                   object=InObject+WHO;
                   char buffer[7]; sprintf(buffer,"%d",WHERE); 
                   line=AtLine+std::string(buffer);
                  }

        virtual const char* which(void) const { return function.data();}
        virtual const char* who(void) const { return object.data();}
        virtual const char* where(void) const { return line.data();}

        virtual void Message(void){ std::cout<<"\n"<<what()<<" "<<which()<<" "<<who()<<" "<<where()<<std::flush;}
       };

// ****************************************************************************

class CANNOT_FIND : public std::runtime_error, public BASIC_ERROR
{ public : CANNOT_FIND(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::runtime_error(std::string("Cannot Find")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           CANNOT_FIND(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::runtime_error(std::string("Cannot Find")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~CANNOT_FIND(void) throw() {;}
};

class CANNOT_READ : public std::runtime_error, public BASIC_ERROR
{ public : CANNOT_READ(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::runtime_error(std::string("Cannot Read")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           CANNOT_READ(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::runtime_error(std::string("Cannot Read")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~CANNOT_READ(void) throw() {;}
};

class CANNOT_WRITE : public std::runtime_error, public BASIC_ERROR
{ public : CANNOT_WRITE(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::runtime_error(std::string("Cannot Write")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           CANNOT_WRITE(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::runtime_error(std::string("Cannot Write")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~CANNOT_WRITE(void) throw() {;}
};

// ****************************************************************************

class EMPTY : public std::length_error, public BASIC_ERROR
{ public : EMPTY(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::length_error(std::string("The object is empty")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           EMPTY(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::length_error(std::string("The object is empty")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~EMPTY(void) throw() {;}
};

class DIFFERENT_DIMENSIONS : public std::length_error, public BASIC_ERROR
{ public : DIFFERENT_DIMENSIONS(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::length_error(std::string("Different dimensions")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           DIFFERENT_DIMENSIONS(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::length_error(std::string("Different dimensions")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~DIFFERENT_DIMENSIONS(void) throw() {;}
};

class DIFFERENT_LENGTHS : public std::length_error, public BASIC_ERROR
{ public : DIFFERENT_LENGTHS(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::length_error(std::string("Different lengths")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           DIFFERENT_LENGTHS(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::length_error(std::string("Different lengths")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~DIFFERENT_LENGTHS(void) throw() {;}
};

class DIFFERENT_SIZES : public std::length_error, public BASIC_ERROR
{ public : DIFFERENT_SIZES(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::length_error(std::string("Different sizes")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           DIFFERENT_SIZES(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::length_error(std::string("Different sizes")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~DIFFERENT_SIZES(void) throw() {;}
};

// ****************************************************************************

class DIVISION_BY_ZERO : public std::invalid_argument, public BASIC_ERROR
{ public : DIVISION_BY_ZERO(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(std::string("Division by zero")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           DIVISION_BY_ZERO(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::invalid_argument(std::string("Division by zero")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~DIVISION_BY_ZERO(void) throw() {;}
};

class EVEN_NUMBER : public std::invalid_argument, public BASIC_ERROR
{ public : EVEN_NUMBER(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(std::string("The number was even")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           EVEN_NUMBER(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::invalid_argument(std::string("The number was even")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~EVEN_NUMBER(void) throw() {;}
};

class NEGATIVE_INTEGER : public std::invalid_argument, public BASIC_ERROR
{ public : NEGATIVE_INTEGER(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(std::string("The number was a negative integer")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           NEGATIVE_INTEGER(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::invalid_argument(std::string("The number was a negative integer")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~NEGATIVE_INTEGER(void) throw() {;}
};

class NEGATIVE_NUMBER : public std::invalid_argument, public BASIC_ERROR
{ public : NEGATIVE_NUMBER(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(std::string("The number was negative")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           NEGATIVE_NUMBER(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::invalid_argument(std::string("The number was negative")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~NEGATIVE_NUMBER(void) throw() {;}
};

class NOT_AN_INTEGER : public std::invalid_argument, public BASIC_ERROR
{ public : NOT_AN_INTEGER(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(std::string("The number was not an integer")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           NOT_AN_INTEGER(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::invalid_argument(std::string("The number was not an integer")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~NOT_AN_INTEGER(void) throw() {;}
};

class ODD_NUMBER : public std::invalid_argument, public BASIC_ERROR
{ public : ODD_NUMBER(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(std::string("The number was odd")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ODD_NUMBER(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::invalid_argument(std::string("The number was odd")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~ODD_NUMBER(void) throw() {;}
};

class ZERO_NUMBER : public std::invalid_argument, public BASIC_ERROR
{ public : ZERO_NUMBER(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(std::string("The number is zero")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ZERO_NUMBER(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::invalid_argument(std::string("The number is zero")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~ZERO_NUMBER(void) throw() {;}
};

class FUNCTION_ERROR : public std::invalid_argument, public BASIC_ERROR
{ public : FUNCTION_ERROR(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::invalid_argument(std::string("Cannot evaluate the function")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           FUNCTION_ERROR(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::invalid_argument(std::string("Cannot evaluate the function")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~FUNCTION_ERROR(void) throw() {;}
};

class NO_SOLUTION : public std::runtime_error, public BASIC_ERROR
{ public : NO_SOLUTION(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::runtime_error(std::string("There is no solution")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           NO_SOLUTION(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::runtime_error(std::string("There is no solution")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~NO_SOLUTION(void) throw() {;}
};

// ****************************************************************************

class UNKNOWN_ERROR : public std::runtime_error, public BASIC_ERROR
{ public : UNKNOWN_ERROR(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
           : std::runtime_error(std::string("Unknown Error")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           UNKNOWN_ERROR(const std::string &WHICH,const std::string &WHO,const int &WHERE)
           : std::runtime_error(std::string("Unknown Error")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
           ~UNKNOWN_ERROR(void) throw() {;}
};

// ****************************************************************************

template <class Type> class EQUAL_VALUES : public std::invalid_argument, public BASIC_ERROR
      { private : Type number;
        public : EQUAL_VALUES(Type NUMBER,const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
                 : std::invalid_argument(std::string("The two numbers are equal")), BASIC_ERROR(WHICH,WHO,WHERE){ number=NUMBER;}
                 EQUAL_VALUES(Type NUMBER,const std::string &WHICH,const std::string &WHO,const int &WHERE)
                 : std::invalid_argument(std::string("The two numbers are equal")), BASIC_ERROR(WHICH,WHO,WHERE){ number=NUMBER;}
                 ~EQUAL_VALUES(void) throw() {;}

                 void ChangeValue(Type NUMBER){ number=NUMBER;}
                 Type Number(void){ return number;}
                 void Message(void);
       };

// **************************************************************************

template <class Type> class NOT_EQUAL_VALUES : public std::invalid_argument, public BASIC_ERROR
      { private : Type number1,number2;
        public : NOT_EQUAL_VALUES(Type NUMBER1,Type NUMBER2,const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
                 : std::invalid_argument(std::string("The two numbers are not equal")), BASIC_ERROR(WHICH,WHO,WHERE){ number1=NUMBER1; number2=NUMBER2;}
                 NOT_EQUAL_VALUES(Type NUMBER1,Type NUMBER2,const std::string &WHICH,const std::string &WHO,const int &WHERE)
                 : std::invalid_argument(std::string("The two numbers are not equal")), BASIC_ERROR(WHICH,WHO,WHERE){ number1=NUMBER1; number2=NUMBER2;}
                 ~NOT_EQUAL_VALUES(void) throw() {;}

                 void ChangeValues(Type NUMBER1,Type NUMBER2){ number1=NUMBER1; number2=NUMBER2;}
                 Type FirstNumber(void){ return number1;}
                 Type SecondNumber(void){ return number2;}
                 void Message(void);
       };

// **************************************************************************

template <class Type> class POLE : public std::invalid_argument, public BASIC_ERROR
      { private : Type number;
        public : POLE(Type NUMBER,const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
                 : std::invalid_argument(std::string("There is a pole")), BASIC_ERROR(WHICH,WHO,WHERE){ number=NUMBER;}
                 POLE(Type NUMBER,const std::string &WHICH,const std::string &WHO,const int &WHERE)
                 : std::invalid_argument(std::string("There is a pole")), BASIC_ERROR(WHICH,WHO,WHERE){ number=NUMBER;}
                 ~POLE(void) throw() {;}

                 void ChangeValue(Type NUMBER){ number=NUMBER;}
                 Type Number(void){ return number;}
                 void Message(void);
       };

// **************************************************************************

template <class Type1,class Type2> class INCOMPATABLE_TYPES : public std::invalid_argument, public BASIC_ERROR
      { public : INCOMPATABLE_TYPES(const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
                 : std::invalid_argument(std::string("Incompatable types")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
                 INCOMPATABLE_TYPES(const std::string &WHICH,const std::string &WHO,const int &WHERE)
                 : std::invalid_argument(std::string("Incompatable types")), BASIC_ERROR(WHICH,WHO,WHERE) {;}
                 ~INCOMPATABLE_TYPES(void) throw() {;}

                 void Message(void);
       };

// **************************************************************************

template <class Type> class DID_NOT_CONVERGE : public std::runtime_error, public BASIC_ERROR
      { private : Type valuereached; int iterations;
        public : DID_NOT_CONVERGE(Type VALUEREACHED,int ITERATIONS,const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
                 : std::runtime_error(std::string("Did not converge")), BASIC_ERROR(WHICH,WHO,WHERE){ valuereached=VALUEREACHED; iterations=ITERATIONS;}
                 DID_NOT_CONVERGE(Type VALUEREACHED,int ITERATIONS,const std::string &WHICH,const std::string &WHO,const int &WHERE)
                 : std::runtime_error(std::string("Did not converge")), BASIC_ERROR(WHICH,WHO,WHERE){ valuereached=VALUEREACHED; iterations=ITERATIONS;}
                 ~DID_NOT_CONVERGE(void) throw() {;}

                 Type ValueReached(void){ return valuereached;}
                 void ChangeValueReached(Type VALUEREACHED){ return valuereached=VALUEREACHED;}
                 int Iterations(void){ return iterations;}
                 void Message(void);
       };

// **************************************************************************

template <class Type1,class Type2> class CONVERGED_TO : public std::runtime_error, public BASIC_ERROR
      { private : Type1 valuereached; Type2 interval;
        public : CONVERGED_TO(Type1 VALUEREACHED,Type2 INTERVAL,const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
                 : std::runtime_error(std::string("Converged to")), BASIC_ERROR(WHICH,WHO,WHERE){ valuereached=VALUEREACHED; interval=INTERVAL;}
                 CONVERGED_TO(Type1 VALUEREACHED,Type2 INTERVAL,const std::string &WHICH,const std::string &WHO,const int &WHERE)
                 : std::runtime_error(std::string("Converged to")), BASIC_ERROR(WHICH,WHO,WHERE){ valuereached=VALUEREACHED; interval=INTERVAL;}
                 ~CONVERGED_TO(void) throw() {;}

                 Type1 ValueReached(void){ return valuereached;}
                 Type2 Interval(void){ return interval;}
                 void Message(void);
       };

// **************************************************************************

template <class Type> class OUT_OF_RANGE : public std::out_of_range, public BASIC_ERROR
         { private : Type value, minvalue, maxvalue;
           public : OUT_OF_RANGE(Type VALUE,Type MINVALUE,Type MAXVALUE,const std::string &WHICH=Unknown,const std::string &WHO=Unknown,const std::string &WHERE=Unknown)
                    : std::out_of_range(std::string("Domain Error")), BASIC_ERROR(WHICH,WHO,WHERE){ value=VALUE; minvalue=MINVALUE; maxvalue=MAXVALUE;}
                    OUT_OF_RANGE(Type VALUE,Type MINVALUE,Type MAXVALUE,const std::string &WHICH,const std::string &WHO,const int &WHERE)
                    : std::out_of_range(std::string("Domain Error")), BASIC_ERROR(WHICH,WHO,WHERE){ value=VALUE; minvalue=MINVALUE; maxvalue=MAXVALUE;}
                    ~OUT_OF_RANGE(void) throw() {;}

                    void ChangeValues(Type VALUE,Type MINVALUE,Type MAXVALUE){ value=VALUE; minvalue=MINVALUE; maxvalue=MAXVALUE;}
                    Type Value(void){ return value;}
                    Type LowerLimit(void){ return minvalue;}
                    Type UpperLimit(void){ return maxvalue;}
                    void Message(void);                    
          };

#endif


