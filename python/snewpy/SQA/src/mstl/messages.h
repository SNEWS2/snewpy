#include<iostream>
#include<string>
#include<iomanip>

#include "mstl.h"

#if !defined(BASIC_MESSAGES)
#define BASIC_MESSAGES

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

class BASIC_MESSAGE;

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

//***********************//
//Message Specializations//
//***********************//

class LOADING;
class WRITING;

class USING_APPROXIMATION;

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

#if !defined(_TRACE_LEVEL)
#define _TRACE_LEVEL 99
#endif

//*********************//
//Trace Specializations//
//*********************//

class TRACE;
class ENTERING;
class LEAVING;

//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************

class BASIC_MESSAGE
      { protected : std::string messagestring, function, object;
                    static std::string NoMessage,UnknownFunction,UnknownObject;

        public : BASIC_MESSAGE(const char *MESSAGESTRING=NULL,const char *FUNCTION=NULL,const char *OBJECT=NULL,bool OUTPUT=true);
                 BASIC_MESSAGE(std::string MESSAGESTRING=NoMessage,std::string FUNCTION=UnknownFunction,std::string OBJECT=UnknownObject,bool OUTPUT=true);
                 BASIC_MESSAGE(const BASIC_MESSAGE &BM){ messagestring=BM.messagestring; function=BM.function; object=BM.object;}
                 std::string MessageString(void){ return messagestring;}
                 std::string Function(void){ return function;}
                 std::string Object(void){ return object;}

                 void Message(void);
       };

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

//***********************//
//Message Specializations//
//***********************//

class LOADING : protected BASIC_MESSAGE
{ public : LOADING(const char *FUNCTION,const char *OBJECT=NULL,bool OUTPUT=true) : BASIC_MESSAGE("Loading a file",FUNCTION,OBJECT,OUTPUT) {;}
           LOADING(std::string FUNCTION=UnknownFunction,std::string OBJECT=UnknownObject,bool OUTPUT=true) : BASIC_MESSAGE(std::string("Loading a file"),FUNCTION,OBJECT,OUTPUT) {;}
};
class WRITING : protected BASIC_MESSAGE
{ public : WRITING(const char *FUNCTION,const char *OBJECT=NULL,bool OUTPUT=true) : BASIC_MESSAGE("Writing a file",FUNCTION,OBJECT,OUTPUT) {;}
           WRITING(std::string FUNCTION=UnknownFunction,std::string OBJECT=UnknownObject,bool OUTPUT=true) : BASIC_MESSAGE(std::string("Writing a file"),FUNCTION,OBJECT,OUTPUT) {;}
};
class USING_APPROXIMATION : protected BASIC_MESSAGE
{ public : USING_APPROXIMATION(const char *FUNCTION,const char *OBJECT=NULL,bool OUTPUT=true) : BASIC_MESSAGE("Using an approximation",FUNCTION,OBJECT,OUTPUT) {;}
           USING_APPROXIMATION(std::string FUNCTION=UnknownFunction,std::string OBJECT=UnknownObject,bool OUTPUT=true) : BASIC_MESSAGE(std::string("Using an approximation"),FUNCTION,OBJECT,OUTPUT) {;}
};

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

//*********************//
//Trace Specializations//
//*********************//

class TRACE : protected BASIC_MESSAGE
{ public : TRACE(const char *TRACESTRING) : BASIC_MESSAGE(TRACESTRING) {;}
           TRACE(std::string TRACESTRING) : BASIC_MESSAGE(TRACESTRING) {;}
};

class ENTERING : protected BASIC_MESSAGE
{ public : ENTERING(const char *FUNCTION,const char *OBJECT=NULL) : BASIC_MESSAGE("Entering",FUNCTION,OBJECT) {;}
           ENTERING(std::string FUNCTION,std::string OBJECT=UnknownObject) : BASIC_MESSAGE(std::string("Entering"),FUNCTION,OBJECT) {;}
};
class LEAVING : protected BASIC_MESSAGE
{ public : LEAVING(const char *FUNCTION,const char *OBJECT=NULL) : BASIC_MESSAGE("\tLeaving",FUNCTION,OBJECT) {;}
           LEAVING(std::string FUNCTION,std::string OBJECT=UnknownObject) : BASIC_MESSAGE(std::string("\tLeaving"),FUNCTION,OBJECT) {;}
};

#endif

