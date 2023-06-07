#include "interpolation data.h"

namespace interpolation{

using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;

using std::string;

using std::cout;
using std::flush;
using std::string;
using std::vector;

// *******************************************************
// *******************************************************
// *******************************************************

string XDATA_SINGLESET::CLASS_NAME("XDATA_SINGLESET");

// *******************************************************

XDATA_SINGLESET::XDATA_SINGLESET(int NX,const double *X)
       { CreateX(NX);
         int i;
         for(i=1;i<=NX;i++){ _pSetX(i,X[i-1]);}
        }

// *******************************************************

double XDATA_SINGLESET::X(int i) const
       { if(i<1 || i>NX()){ throw OUT_OF_RANGE<int>(i,1,NX(),string("X"),CLASS_NAME);}
         return _pX(i);
        }

void XDATA_SINGLESET::SetX(int i,double X) const
       { if(i<1 || i>NX()){ throw OUT_OF_RANGE<int>(i,1,NX(),string("SetX"),CLASS_NAME);}
         return _pSetX(i,X);
        }

// *******************************************************

void XDATA_SINGLESET::AddX(double X)
     { vector<double> newx(_pX()); 
       newx.push_back(X);
       _pSetX(newx);
      }

void XDATA_SINGLESET::AddX(const vector<double> &X)
     { vector<double> newx(_pX()); 
       for(int i=1;i<=static_cast<int>(X.size());i++){ newx.push_back(X[i-1]);}
       _pSetX(newx);
      }

// *******************************************************

// shift, rescale the data      
void XDATA_SINGLESET::ShiftX(const double &S)
     { if(Empty()==true){ throw EMPTY("ShiftX",CLASS_NAME,77);}

       vector<double> XT(NX());
       int i;
       for(i=1;i<=NX();i++){ XT[i-1]=_pX(i)+S;}
       _pSetX(XT);
      }

void XDATA_SINGLESET::RescaleX(const double &S)
     { if(Empty()==true){ throw EMPTY("RescaleX",CLASS_NAME,85);}

       vector<double> XT(NX());
       int i;
       for(i=1;i<=NX();i++){ XT[i-1]=_pX(i)*S;}
       _pSetX(XT);
      }

void XDATA_SINGLESET::AbsX(void)
          { int i;
            for(i=1;i<=NX();i++){ _pSetX(i,Sign(_pX(i))*_pX(i));} 
           }

// *******************************************************
// *******************************************************
// *******************************************************

string YDATA_SINGLESET<1>::CLASS_NAME("YDATA_SINGLESET<1>");

// *******************************************************

YDATA_SINGLESET<1>::YDATA_SINGLESET(int NY,const double *Y)
       { CreateY(NY);
         int i;
         for(i=1;i<=NY;i++){ _pSetY(i,Y[i-1]);}
        }

// *******************************************************

double YDATA_SINGLESET<1>::Y(int i) const
       { if(i<1 || i>NY()){ throw OUT_OF_RANGE<int>(i,1,NY(),string("Y"),CLASS_NAME);}
         return _pY(i);
        }

void YDATA_SINGLESET<1>::SetY(int i,double Y) const 
       { if(i<1 || i>NY()){ throw OUT_OF_RANGE<int>(i,1,NY(),string("SetY"),CLASS_NAME);}
         return _pSetY(i,Y);
        }

// *******************************************************

void YDATA_SINGLESET<1>::AddY(double Y)
     { vector<double> newy(_pY()); 
       newy.push_back(Y);
       _pSetY(newy);
      }

void YDATA_SINGLESET<1>::AddY(const vector<double> &Y)
     { vector<double> newy(_pY()); 
       for(int i=1;i<=static_cast<int>(Y.size());i++){ newy.push_back(Y[i-1]);}
       _pSetY(newy);
      }

// *******************************************************

// shift, rescale the data      
void YDATA_SINGLESET<1>::ShiftY(const double &S)
     { if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME,77);}

       vector<double> YT(NY());
       int i;
       for(i=1;i<=NY();i++){ YT[i-1]=_pY(i)+S;}
       _pSetY(YT);
      }

void YDATA_SINGLESET<1>::RescaleY(const double &S)
     { if(Empty()==true){ throw EMPTY("RescaleY",CLASS_NAME,85);}

       vector<double> YT(NY());
       int i;
       for(i=1;i<=NY();i++){ YT[i-1]=_pY(i)*S;}
       _pSetY(YT);
      }

void YDATA_SINGLESET<1>::AbsY(void)
          { int i;
            for(i=1;i<=NY();i++){ _pSetY(i,Sign(_pY(i))*_pY(i));} 
           }

// *******************************************************
// *******************************************************
// *******************************************************

void DELTAX_SINGLESET::XDifference(void) const
     { CreateDeltaX(NX()-1);
       int a;
       for(a=1;a<=NX()-1;a++){ _pSetDeltaX(a,_pX(a+1)-_pX(a));}
       SetXDifferenced(true);
      }

// *******************************************************
// *******************************************************
// *******************************************************

void DELTAY_SINGLESET::YDifference(void) const
     { CreateDeltaY(NY()-1);
       int a;
       for(a=1;a<=NY()-1;a++){ _pSetDeltaY(a,_pY(a+1)-_pY(a));}
       SetYDifferenced(true);
      }

// *******************************************************
// *******************************************************
// *******************************************************

void XLIMITS_SINGLESET::XLimit(void) const
       { if(Empty()==true){ throw EMPTY("XLimit");}
          double xmin=_pX(1);
          double xmax=_pX(1);
          for(int i=2;i<=NX();i++)
             { if(_pX(i)<xmin){ xmin=_pX(i);}
                if(_pX(i)>xmax){ xmax=_pX(i);}
              }
         SetXMin(xmin);
         SetXMax(xmax);
         SetXLimited(true);
        }

double XLIMITS_SINGLESET::XMin(void) const { if(XLimited()==false){ XLimit();} return _pXMin();}

double XLIMITS_SINGLESET::XMax(void) const { if(XLimited()==false){ XLimit();} return _pXMax();}

// *******************************************************

void YLIMITS_SINGLESET::YLimit(void) const
       { if(Empty()==true){ throw EMPTY("YLimit");}
          double ymin=_pY(1);
         double ymax=_pY(1);
         for(int i=2;i<=NY();i++)
            { if(_pY(i)<ymin){ ymin=_pY(i);}
              if(_pY(i)>ymax){ ymax=_pY(i);}
             }
         SetYMin(ymin);
         SetYMax(ymax);
         SetYLimited(true);
        }

double YLIMITS_SINGLESET::YMin(void) const { if(YLimited()==false){ YLimit();} return _pYMin();}

double YLIMITS_SINGLESET::YMax(void) const { if(YLimited()==false){ YLimit();} return _pYMax();}

// *******************************************************

void XYLIMITS_SINGLEXSET_SINGLEYSET::XYLimit(void) const
          { double xmin=X(1), xmax=X(1), yatxmin=Y(1), yatxmax=Y(1);
             double ymin=Y(1), ymax=Y(1), xatymin=X(1), xatymax=X(1);

             for(int i=2;i<=NX();i++)
                   { if(X(i)<xmin){ xmin=X(i); yatxmin=Y(i);}
                      if(X(i)>xmax){ xmax=X(i); yatxmax=Y(i);}
                      if(Y(i)<ymin){ ymin=Y(i); xatymin=X(i);}
                      if(Y(i)>ymax){ ymax=Y(i); xatymax=X(i);}
                     }
         SetYAtXMin(yatxmin);
         SetYAtXMax(yatxmax);
         SetXAtYMin(xatymin);
         SetXAtYMax(xatymax);

         SetXYLimited(true);
        }

double XYLIMITS_SINGLEXSET_SINGLEYSET::YAtXMin(void) const { if(XYLimited()==false){ XYLimit();} return _pYAtXMin();}
double XYLIMITS_SINGLEXSET_SINGLEYSET::YAtXMax(void) const { if(XYLimited()==false){ XYLimit();} return _pYAtXMax();}

double XYLIMITS_SINGLEXSET_SINGLEYSET::XAtYMin(void) const { if(XYLimited()==false){ XYLimit();} return _pXAtYMin();}
double XYLIMITS_SINGLEXSET_SINGLEYSET::XAtYMax(void) const { if(XYLimited()==false){ XYLimit();} return _pXAtYMax();}

// *******************************************************
// *******************************************************
// *******************************************************
// *******************************************************

string GRADEBASE_SINGLEXSET_SINGLEYSET::CLASS_NAME("GRADEBASE_SINGLEXSET_SINGLEYSET");

// **************

double GRADEBASE_SINGLEXSET_SINGLEYSET::G(int i) const
       { if(i<1 || i>NG()){ throw OUT_OF_RANGE<int>(i,1,NG(),string("G"),CLASS_NAME);}
         return _pG(i);
        }

void GRADEBASE_SINGLEXSET_SINGLEYSET::SetG(int i,double G) const
       { if(i<1 || i>NG()){ throw OUT_OF_RANGE<int>(i,1,NG(),string("SetG"),CLASS_NAME);}
         return _pSetG(i,G);
        }

// *******************************************************************

void LINEGRADE_SINGLEXSET_SINGLEYSET::Grade(void) const
     { CreateG(NX());
       int i;
               for(i=1;i<=NG()-1;i++){ _pSetG( i,_pDeltaY(i)/_pDeltaX(i) );}
       _pSetG( NG(),_pG(NG()-1) );
       
       SetGraded(true);       
      }

// *******************************************************

void LOCALGRADE_SINGLEXSET_SINGLEYSET::Grade(void) const
     { CreateG(NX());
       if(NG()<=5){ throw TOO_FEW_POINTS("Grade","LOCALGRADE_SINGLEXSET_SINGLEYSET");}

       vector<double> XX(5), YY(5), D;

       int i,xoffset;
       for(i=1;i<=NG();i++)
             { if(i==1){ xoffset=0;} 
                else{ if(i==2){ xoffset=-1;}  
                            else{ if(i==NG()-1){ xoffset=-3;} 
                                        else{ if(i==NG()){ xoffset=-4;} 
                                                    else{ xoffset=-2;} 
                                                   } 
                                      } 
                            }
                         
                for(int a=0;a<=4;a++){ XX[a]=_pX(i+a+xoffset)-_pX(i); YY[a]=_pY(i+a+xoffset)-_pY(i);}
               D=FiniteDifference1D(0.,XX,YY);
               _pSetG(i,D[1]);
              }

       SetGraded(true);
      }

// *******************************************************
// *******************************************************
// *******************************************************

void SPLINE_SINGLEXSET_SINGLEYSET::DestroyA(void) const
      { if(Fitted()==false){ return;}
        int i;
        for(i=1;i<=NA();i++){ a[i-1].clear();} 
        SetFitted(false);
       }

// *************

double SPLINE_SINGLEXSET_SINGLEYSET::Interpolate(int i,double T) const
       { double I=A(i,NParameters(i)-1);
         for(int j=NParameters(i)-2;j>=0;j--){ (I*=T)+=A(i,j);}
         return I;
        }

double SPLINE_SINGLEXSET_SINGLEYSET::Derivative(int i,double T) const
       { double D=A(i,NParameters(i)-1)*(NParameters(i)-1.);
         for(int j=NParameters(i)-2;j>=1;j--){ (D*=T)+=A(i,j)*j;}
         return D;
        }

double SPLINE_SINGLEXSET_SINGLEYSET::SecondDerivative(int i,double T) const
       { double SD=A(i,NParameters(i)-1)*(NParameters(i)-1.)*(NParameters(i)-2.);
         for(int j=NParameters(i)-2;j>=2;j--){ (SD*=T)+=A(i,j)*j*(j-1);}
         return SD;
        }

double SPLINE_SINGLEXSET_SINGLEYSET::Integral(int i,double T0,double T1) const
       { double I1,I0;
         I0=I1=A(i,NParameters(i)-1)/(NParameters(i)-1);
         for(int j=NParameters(i)-2;j>=0;j--){ (I0*=T0)+=A(i,j)/(j+1); (I1*=T1)+=A(i,j)/(j+1);}
         I0*=T0; I1*=T1;
         return I1-I0;
        }

// *******************************************************************
// *******************************************************************
// *******************************************************************

void SORT_SINGLEXSET_SINGLEYSET::Sort(void)
     { double tmpx, tmpy;
       for(int j=2;j<=NX();j++)
          { tmpx=_pX(j); tmpy=_pY(j);
            int i=j-1;
            while(i>=1 && _pX(i)>tmpx){ _pSetX(i+1,_pX(i)); _pSetY(i+1,_pY(i)); i--;};
            _pSetX(i+1,tmpx); _pSetY(i+1,tmpy);       
           }
       SetSorted(true); 
      }

// *******************************************************************
// *******************************************************************
// *******************************************************************

void OPENWRITE_SINGLEXSET_SINGLEYSET::Open(string filename)
       { ifstream fin(filename.c_str());

         if(fin){ int N, offset;
                  double d; 
                  vector<double> data;             
             
                  LOADING("Open",filename);
                  while(fin>>d){ data.push_back(d);};
                  fin.close();
             
                  if(Odd(data.size())==true){ N=static_cast<int>(data[0]); offset=1;} else{ N=static_cast<int>(data.size()/2); offset=0;}
                  CreateX(N);
                  CreateY(N);
  		  for(int a=0;a<=N-1;a++){ _pSetX(a+1,data[2*a+offset]); _pSetY(a+1,data[2*a+offset+1]);}
                 }
         else{ throw CANNOT_FIND("Open",filename);}
        }

void OPENWRITE_SINGLEXSET_SINGLEYSET::Open(string filename,char ignore)
       { ifstream fin(filename.c_str());

         if(fin){ int N, offset;
                  double d; 
                  vector<double> data;             
                  string line;
             
                  LOADING("Open",filename);    
                  while(fin.peek()==ignore){ getline(fin,line);};
                  while(fin>>d){ data.push_back(d);};
                  fin.close();
             
                  if(Odd(data.size())==true){ N=static_cast<int>(data[0]); offset=1;} else{ N=static_cast<int>(data.size()/2); offset=0;}
                  CreateX(N);
                  CreateY(N);
  		  for(int a=0;a<=N-1;a++){ _pSetX(a+1,data[2*a+offset]); _pSetY(a+1,data[2*a+offset+1]);}
                 }
         else{ throw CANNOT_FIND("Open",filename);}
        }

void OPENWRITE_SINGLEXSET_SINGLEYSET::Write(string filename)
     { ofstream fout(filename.c_str());
       if(fout){ WRITING("Write",filename); 
                     fout<<NX();
                     for(int i=1;i<=NX();i++){ fout<<_pX(i)<<"\t"<<_pY(i)<<"\n";}
                     fout.flush();
                     fout.close();
                   }
       else{ throw CANNOT_WRITE("Write",filename);}
      } 

// output the x and y's in point format
ostream& OPENWRITE_SINGLEXSET_SINGLEYSET::operator<<(ostream &os)
         { os<<NX()<<flush; os.precision(15);
           for(int a=1;a<=NX();a++){ os<<"\n"<<_pX(a)<<"\t"<<_pY(a)<<flush;}
           return os;
          }

istream& OPENWRITE_SINGLEXSET_SINGLEYSET::operator>>(istream &is)
         { int N, offset;
           double d; 
           vector<double> data;             
             
           while(is>>d){ data.push_back(d);};
             
           if(Odd(data.size())==true){ N=static_cast<int>(data[0]); offset=1;} else{ N=static_cast<int>(data.size()/2); offset=0;}
           CreateX(N);
           CreateY(N);

           for(int a=0;a<=N-1;a++){ _pSetX(a+1,data[2*a+offset]); _pSetY(a+1,data[2*a+offset+1]);}

           return is;
          }

// *******************************************************************
// *******************************************************************
// *******************************************************************

int XINTERVAL_SINGLESET::XInterval(double X) const
      { if(X<XMin() || X>XMax()){ throw OUT_OF_RANGE<double>(X,XMin(),XMax(),string("XInterval"));}

        if(interval<1 || interval>NX()-1){ interval=(NX()-1)/2;}

        if(interval!=NX()-1 && X>=_pX(interval) && _pX(interval+1)>X){ return interval;}
        if(interval==NX()-1 && X>=_pX(interval) && _pX(interval+1)>=X){ return interval;}

        if(interval>=2){ if(X>=_pX(interval-1) && _pX(interval)>X){ --interval; return interval;} }
        if(interval<=NX()-3){ if(X>=_pX(interval+1) && _pX(interval+2)>X){ ++interval; return interval;} }
        if(interval==NX()-2){ if(X>=_pX(interval+1) && _pX(interval+2)>=X){ ++interval; return interval;} }
  
        if(Equality(X,_pX(1))==true){ return interval=1;} if(Equality(X,_pX(NX()))==true){ return interval=NX()-1;}

        int lower=1, upper=NX()-1;
        interval=( upper+lower + (upper-lower)%2 )/2;

        while(lower!=upper && (X<_pX(interval) || X>=_pX(interval+1)) )
             { if(X>=_pX(interval+1)){ lower=interval+1;} else{ upper=interval-1;}
               interval=( upper+lower + (upper-lower)%2 )/2;
              }

        return interval;
       }

// *******************************************************************

int YINTERVAL_SINGLESET::YInterval(double Y) const
    { if(Y<YMin() || Y>YMax()){ throw OUT_OF_RANGE<double>(Y,YMin(),YMax(),string("YInterval"));}

      if(interval<1 || interval>NY()-1){ interval=(NY()-1)/2;}

      if(interval!=NY()-1 && Y>=_pY(interval) && _pY(interval+1)>Y){ return interval;}
      if(interval==NY()-1 && Y>=_pY(interval) && _pY(interval+1)>=Y){ return interval;}

      if(interval>=2){ if(Y>=_pY(interval-1) && _pY(interval)>Y){ --interval; return interval;} }
      if(interval<=NY()-3){ if(Y>=_pY(interval+1) && _pY(interval+2)>Y){ ++interval; return interval;} }
      if(interval==NY()-2){ if(Y>=_pY(interval+1) && _pY(interval+2)>=Y){ ++interval; return interval;} }

      if(Equality(Y,_pY(1))==true){ return interval=1;} if(Equality(Y,_pY(NY()))==true){ return interval=NY()-1;}

      int lower=1, upper=NY()-1;
      interval=( upper+lower + (upper-lower)%2 )/2;

      while(lower!=upper && (Y<_pY(interval) || Y>=_pY(interval+1)) )
           { if(Y>=_pY(interval+1)){ lower=interval+1;} else{ upper=interval-1;}
        interval=( upper+lower + (upper-lower)%2 )/2;
       }

      return interval; 
     }

// *******************************************************************
// *******************************************************************
// *******************************************************************

void DISTINGUISH_SINGLEYSET::Distinguish(void)
     { if(Sorted()==false){ Sort();}
       cout.precision(7); 

       int i; 
       vector<double> newx, newy; 
       newx.reserve(N()); 
       newy.reserve(N());

       for(i=1;i<=N()-1;i++)
          { if(Equality(_pX(i),_pX(i+1))==false){ newx.push_back(_pX(i)); newy.push_back(_pY(i));}
            else{ cout<<"\nReplacing the point "<<i<<" with x= "<<_pX(i)<<" y= "<<_pY(i);
                  int j=1; double mean=_pY(i); bool end=false;
                  while(end==false){ cout<<"\nand the point "<<i+j<<" with x= "<<_pX(i+j)<<" y= "<<_pY(i+j);
                           mean=( mean*j + _pY(i+j) )/(1.+j);
                           if(i+j>=N()){ end=true;} else{ if(Equality(_pX(i),_pX(i+j+1))==false){ end=true;} else{ j++;} }
                          };
                  cout<<"\nwith a point at x= "<<_pX(i)<<" and y= "<<mean; cout.flush();
                  newx.push_back(_pX(i)); newy.push_back(mean);
                  i+=j;
                 }
           }
       if(i==N()){ newx.push_back(_pX(N())); newy.push_back(_pY(N()));}
     
       _pSetX(newx);
       _pSetY(newy); 
      }

// *******************************************************
// *******************************************************
// *******************************************************

OPERATOR& OPERATOR::operator+=(const double &A)
            { int i;
              for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)+A );} 
              return (*this);
             }
OPERATOR& OPERATOR::operator-=(const double &A)
             { int i;
               for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)-A );} 
               return (*this);
             }
OPERATOR& OPERATOR::operator*=(const double &A)
             { int i;
               for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)*A );} 
               return (*this);
              }
OPERATOR& OPERATOR::operator/=(const double &A)
       { if(Equality(A,0.)==true){ throw DIVISION_BY_ZERO("/=double");}
         int i;
         for(i=1;i<=NY();i++){ _pSetY(i,_pY(i)/A);}
         return (*this);
        }

// *******************************************************************

OPERATOR& OPERATOR::operator=(const NULLARYFUNCTOR<double> &NF)
            { int i;
              for(i=1;i<=NY();i++){ _pSetY( i,NF() );} 
              return (*this);
             }
OPERATOR& OPERATOR::operator+=(const NULLARYFUNCTOR<double> &NF)
            { int i;
              for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)+NF() );}
              return (*this);
             }
OPERATOR& OPERATOR::operator-=(const NULLARYFUNCTOR<double> &NF)
             { int i;
               for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)-NF() );} 
               return (*this);
             }
OPERATOR& OPERATOR::operator*=(const NULLARYFUNCTOR<double> &NF)
             { int i;
               for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)*NF() );} 
               return (*this);
              }
OPERATOR& OPERATOR::operator/=(const NULLARYFUNCTOR<double> &NF)
       { double D;
         int i;
         for(i=1;i<=NY();i++)
            { D=NF();
               if(Equality(D,0.)==true){ throw DIVISION_BY_ZERO("/=NULLARYFUNCTOR");}
              _pSetY(i,_pY(i)/D);
             }
         return (*this);
        }

// *******************************************************************

OPERATOR& OPERATOR::operator=(const UNARYFUNCTOR<double> &UF)
             { int i;
               for(i=1;i<=NY();i++){ _pSetY( i,UF(_pX(i)) );} 
               return (*this);
              }
OPERATOR& OPERATOR::operator+=(const UNARYFUNCTOR<double> &UF)
             { int i;
               for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)+UF(_pX(i)) );} 
               return (*this);
             }
OPERATOR& OPERATOR::operator-=(const UNARYFUNCTOR<double> &UF)
             { int i;
               for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)-UF(_pX(i)) );} 
               return (*this);
              }
OPERATOR& OPERATOR::operator*=(const UNARYFUNCTOR<double> &UF)
            { int i;
              for(i=1;i<=NY();i++){ _pSetY( i,_pY(i)*UF(_pX(i)) );} 
              return (*this);
             }
OPERATOR& OPERATOR::operator/=(const UNARYFUNCTOR<double> &UF)
        { double D;
          int i;
          for(i=1;i<=NY();i++)
                 { D=UF(_pX(i));
                     if(Equality(D,0.)==true){ throw DIVISION_BY_ZERO("/=UNARYFUNCTOR");}
                      _pSetY(i,_pY(i)/D);
                   }
          return (*this);
        }

// *******************************************************************

OPERATOR& OPERATOR::operator+=(const BINARYFUNCTOR<double> &BF)
            { int i;
              for(i=1;i<=NY();i++){ _pSetY(i,_pY(i)+BF(_pX(i),_pY(i)));} 
              return (*this);
             }
OPERATOR& OPERATOR::operator-=(const BINARYFUNCTOR<double> &BF)
            { int i;
              for(i=1;i<=NY();i++){ _pSetY(i,_pY(i)-BF(_pX(i),_pY(i)));} 
              return (*this);
             }
OPERATOR& OPERATOR::operator*=(const BINARYFUNCTOR<double> &BF)
            { int i;
              for(i=1;i<=NY();i++){ _pSetY(i,_pY(i)*BF(_pX(i),_pY(i)));} 
              return (*this);
             }
OPERATOR& OPERATOR::operator/=(const BINARYFUNCTOR<double> &BF)
       { double D;
         int i;
         for(i=1;i<=NY();i++)
                 { D=BF(_pX(i),_pY(i));
                    if(Equality(D,0.)==true){ throw DIVISION_BY_ZERO("/=BINARYFUNCTOR");}
                    _pSetY(i,_pY(i)/D);
                   }
         return (*this);
        }

} //end of namespace
