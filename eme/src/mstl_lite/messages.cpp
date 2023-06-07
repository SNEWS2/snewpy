#include "messages.h"

using std::cout;
using std::string;
using std::flush;
using std::endl;

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

string BASIC_MESSAGE::NoMessage("No Message");
string BASIC_MESSAGE::UnknownFunction("Unkown Function");
string BASIC_MESSAGE::UnknownObject("Unkown Object");

BASIC_MESSAGE::BASIC_MESSAGE(const char *MESSAGESTRING,const char *FUNCTION,const char *OBJECT,bool OUTPUT)
               { if(MESSAGESTRING!=NULL){ messagestring=string(MESSAGESTRING);}
                 if(FUNCTION!=NULL){ function=string(FUNCTION);}
                 if(OBJECT!=NULL){ object=string(OBJECT);}
                 if(OUTPUT==true){ Message();}
                }

BASIC_MESSAGE::BASIC_MESSAGE(string MESSAGESTRING,string FUNCTION,string OBJECT,bool OUTPUT)
               { messagestring=string(MESSAGESTRING);
                 function=string(FUNCTION);
                 object=string(OBJECT);
                 if(OUTPUT==true){ Message();}
                }

void BASIC_MESSAGE::Message(void)
     { cout<<"\n"<<MessageString();
       if(function.empty()==false){ cout<<" Thrown in the function "<<Function();}
       if(object.empty()==false){ cout<<" by the object "<<Object();}
       cout<<flush;
      }

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

