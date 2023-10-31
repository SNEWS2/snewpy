
#include "miscellaneous functions.h"

//*****************************************************************************************
//*****************************************************************************************
//*****************************************************************************************

void hold(int s){ time_t t0=time(NULL),t1; do{ t1=time(NULL);} while(t1<s+t0);}

void DoNothing(void){;}

bool Odd(int i){ return i&1;}
bool Even(int i){ return !Odd(i);}

void Zero(void){;}

void One(void){;}

void Two(void){;}
