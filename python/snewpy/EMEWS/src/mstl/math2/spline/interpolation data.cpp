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

XDATA_SINGLESET::XDATA_SINGLESET(size_t NX,const double *X)
       { CreateX(NX);
         for(int i=1;i<=static_cast<int>(NX);i++){ _pSetX(i,X[i-1]);}
        }

// *******************************************************

double XDATA_SINGLESET::X(int i) const
       { if(i<1 || i>static_cast<int>(NX())){ throw OUT_OF_RANGE<int>(i,1,NX(),string("X"),CLASS_NAME);}
         return _pX(i);
        }

void XDATA_SINGLESET::SetX(int i,double X) const
       { if(i<1 || i>static_cast<int>(NX())){ throw OUT_OF_RANGE<int>(i,1,NX(),string("SetX"),CLASS_NAME);}
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
       for(int i=1;i<=static_cast<int>(NX());i++){ XT[i-1]=_pX(i)+S;}
       _pSetX(XT);
      }

void XDATA_SINGLESET::RescaleX(const double &S)
     { if(Empty()==true){ throw EMPTY("RescaleX",CLASS_NAME,85);}

       vector<double> XT(NX());
       for(int i=1;i<=static_cast<int>(NX());i++){ XT[i-1]=_pX(i)*S;}
       _pSetX(XT);
      }

void XDATA_SINGLESET::AbsX(void)
          { for(int i=1;i<=static_cast<int>(NX());i++){ _pSetX(i,Sign(_pX(i))*_pX(i));} 
           }

// *******************************************************
// *******************************************************
// *******************************************************

string XDATA_MULTIPLESETS::CLASS_NAME("XDATA_MULTIPLESETS");

// *******************************************************

void XDATA_MULTIPLESETS::CreateX(size_t NXSets,std::vector<size_t> NX) const
          { CreateX(NXSets);
            for(int i=1;i<=static_cast<int>(NXSets);i++){ x[i-1]=std::vector<double>(NX[i-1]);} 
           }

// *******************************************************

void XDATA_MULTIPLESETS::DestroyX(void) const
          { for(int i=1;i<=static_cast<int>(NXSets());i++){ x[i-1].clear();} 
           }

// *******************************************************

XDATA_MULTIPLESETS::XDATA_MULTIPLESETS(size_t NXSets,size_t NX,const double **X)
       { CreateX(NXSets,NX);
          for(int i=1;i<=static_cast<int>(NXSets);i++){ for(int j=1;j<=static_cast<int>(NX);j++){ _pSetX(i,j,X[i-1][j-1]);} }
        }

// *******************************************************

vector<size_t> XDATA_MULTIPLESETS::NX(void) const
      { vector<size_t> n(NXSets());
        for(int i=1;i<=static_cast<int>(NXSets());i++){ n[i-1]=x[i-1].size();}
        return n;
       }

// *******************************************************

double XDATA_MULTIPLESETS::X(int i,int j) const
       { if(i<1 || i>static_cast<int>(NXSets())){ throw OUT_OF_RANGE<int>(i,1,NXSets(),string("X"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NX(i))){ throw OUT_OF_RANGE<int>(j,1,NX(i),string("X"),CLASS_NAME);}
         return _pX(i,j);
        }

vector<double> XDATA_MULTIPLESETS::X(int i) const
       { if(i<1 || i>static_cast<int>(NXSets())){ throw OUT_OF_RANGE<int>(i,1,NXSets(),string("X"),CLASS_NAME);}
         return _pX(i);
        }

vector<double> XDATA_MULTIPLESETS::X(vector<int> i) const
       { vector<double> xx(NXSets());
         for(int j=1;j<=static_cast<int>(NXSets());j++)
            { if(i[j-1]<1 || i[j-1]>static_cast<int>(NX(j))){ throw OUT_OF_RANGE<int>(i[j-1],1,NX(j),string("X"),CLASS_NAME);}
              xx[j-1]=_pX(j,i[j-1]);
             } 
         return xx;
        }         

void XDATA_MULTIPLESETS::SetX(int i,int j,double X) const
       { if(i<1 || i>static_cast<int>(NXSets())){ throw OUT_OF_RANGE<int>(i,1,NXSets(),string("SetX"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NX(i))){ throw OUT_OF_RANGE<int>(j,1,NX(i),string("SetX"),CLASS_NAME);}
         return _pSetX(i,j,X);
        }

void XDATA_MULTIPLESETS::SetX(int i,const vector<double> &X) const
       { if(i<1 || i>static_cast<int>(NXSets())){ throw OUT_OF_RANGE<int>(i,1,NXSets(),string("SetX"),CLASS_NAME);}
         return _pSetX(i,X);
        }

// *******************************************************

void XDATA_MULTIPLESETS::AddX(int i,const vector<double> &X)
     { vector<double> newx(_pX(i)); 
       for(int j=1;j<=static_cast<int>(X.size());j++){ newx.push_back(X[j-1]);}
       _pSetX(i,newx);
      }

void XDATA_MULTIPLESETS::AddX(const vector<vector<double> > &X)
     { if(X.size()!=NXSets()){ throw DIFFERENT_LENGTHS("AddX",CLASS_NAME,163);}
        for(int i=1;i<=static_cast<int>(NXSets());i++)
             { vector<double> newx(_pX(i));            
                for(int j=1;j<=static_cast<int>(X[i].size());j++){ newx.push_back(X[i-1][j-1]);}
                _pSetX(i,newx);
              }
      }

// *******************************************************

// shift, rescale the data      
void XDATA_MULTIPLESETS::ShiftX(int i,const double &S)
     { if(Empty()==true){ throw EMPTY("ShiftX",CLASS_NAME,176);}
       if(Empty(i)==true){ throw EMPTY("ShiftX",CLASS_NAME,177);}

       vector<double> XT(NX(i));
       for(int j=1;j<=static_cast<int>(NX(i));j++){ XT[j-1]=_pX(i,j)+S;}
       _pSetX(i,XT);
      }

void XDATA_MULTIPLESETS::ShiftX(const double &S)
     { if(Empty()==true){ throw EMPTY("ShiftX",CLASS_NAME,167);}
       for(int i=1;i<=static_cast<int>(NXSets());i++)
          { if(Empty(i)==true){ throw EMPTY("ShiftX",CLASS_NAME,169);}
            vector<double> XT(NX(i));
            for(int j=1;j<=static_cast<int>(NX(i));j++){ XT[j-1]=_pX(i,j)+S;}
            _pSetX(i,XT);
          }
      }

void XDATA_MULTIPLESETS::RescaleX(int i,const double &S)
     { if(Empty()==true){ throw EMPTY("RescaleX",CLASS_NAME,177);}
       if(Empty(i)==true){ throw EMPTY("RescaleX",CLASS_NAME,178);}

       vector<double> XT(NX(i));
       for(int j=1;j<=static_cast<int>(NX(i));j++){ XT[j-1]=_pX(i,j)*S;}
       _pSetX(i,XT);
      }

void XDATA_MULTIPLESETS::RescaleX(const double &S)
     { if(Empty()==true){ throw EMPTY("RescaleX",CLASS_NAME,186);}
       for(int i=1;i<=static_cast<int>(NXSets());i++)
          { if(Empty(i)==true){ throw EMPTY("RescaleX",CLASS_NAME,188);}
            vector<double> XT(NX(i));
            for(int j=1;j<=static_cast<int>(NX(i));j++){ XT[j-1]=_pX(i,j)*S;}
            _pSetX(i,XT);
           }
      }

void XDATA_MULTIPLESETS::AbsX(int i)
          { int j;
            for(j=1;j<=static_cast<int>(NX(i));j++){ _pSetX(i,j,Sign(_pX(i,j))*_pX(i,j));} 
           }

void XDATA_MULTIPLESETS::AbsX(void)
          { for(int i=1;i<=static_cast<int>(NXSets());i++)
                   { for(int j=1;j<=static_cast<int>(NX(i));j++){ _pSetX(i,j,Sign(_pX(i,j))*_pX(i,j));} 
                     } 
            }

// *******************************************************
// *******************************************************
// *******************************************************

string YDATA_SINGLESET<1>::CLASS_NAME("YDATA_SINGLESET<1>");

// *******************************************************

YDATA_SINGLESET<1>::YDATA_SINGLESET(size_t NY,const double *Y)
       { CreateY(NY);
         for(int i=1;i<=static_cast<int>(NY);i++){ _pSetY(i,Y[i-1]);}
        }

// *******************************************************

double YDATA_SINGLESET<1>::Y(int i) const
       { if(i<1 || i>static_cast<int>(NY())){ throw OUT_OF_RANGE<int>(i,1,NY(),string("Y"),CLASS_NAME);}
         return _pY(i);
        }

void YDATA_SINGLESET<1>::SetY(int i,double Y) const 
       { if(i<1 || i>static_cast<int>(NY())){ throw OUT_OF_RANGE<int>(i,1,NY(),string("SetY"),CLASS_NAME);}
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
       for(int i=1;i<=static_cast<int>(NY());i++){ YT[i-1]=_pY(i)+S;}
       _pSetY(YT);
      }

void YDATA_SINGLESET<1>::RescaleY(const double &S)
     { if(Empty()==true){ throw EMPTY("RescaleY",CLASS_NAME,85);}

       vector<double> YT(NY());
       for(int i=1;i<=static_cast<int>(NY());i++){ YT[i-1]=_pY(i)*S;}
       _pSetY(YT);
      }

void YDATA_SINGLESET<1>::AbsY(void)
          { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY(i,Sign(_pY(i))*_pY(i));} 
           }

// *******************************************************
// *******************************************************
// *******************************************************

string YDATA_MULTIPLESETS::CLASS_NAME("YDATA_MULTIPLESETS");

// *******************************************************

void YDATA_MULTIPLESETS::CreateY(size_t NYSets,std::vector<size_t> NY) const
          { CreateY(NYSets);
            for(int i=1;i<=static_cast<int>(NYSets);i++){ y[i-1]=std::vector<double>(NY[i-1]);} 
           }

// *******************************************************

void YDATA_MULTIPLESETS::DestroyY(void) const
     { for(int i=1;i<=static_cast<int>(NYSets());i++){ y[i-1].clear();} 
      }

// *******************************************************

YDATA_MULTIPLESETS::YDATA_MULTIPLESETS(size_t NYSets,size_t NY,const double **Y)
       { CreateY(NYSets,NY);
         for(int i=1;i<=static_cast<int>(NYSets);i++)
               { for(int j=1;j<=static_cast<int>(NY);j++){ _pSetY(i,j,Y[i-1][j-1]);} 
                }
        }

// *******************************************************

vector<size_t> YDATA_MULTIPLESETS::NY(void) const
      { vector<size_t> n(NYSets());
        for(int i=1;i<=static_cast<int>(NYSets());i++){ n[i-1]=y[i-1].size();}
        return n;
       }

// *******************************************************

double YDATA_MULTIPLESETS::Y(int i,int j) const
       { if(i<1 || i>static_cast<int>(NYSets())){ throw OUT_OF_RANGE<int>(i,1,NYSets(),string("Y"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NY(i))){ throw OUT_OF_RANGE<int>(j,1,NY(i),string("Y"),CLASS_NAME);}
         return _pY(i,j);
        }

vector<double> YDATA_MULTIPLESETS::Y(int i) const
       { if(i<1 || i>static_cast<int>(NYSets())){ throw OUT_OF_RANGE<int>(i,1,NYSets(),string("Y"),CLASS_NAME);}
         return _pY(i);
        }

vector<double> YDATA_MULTIPLESETS::Y(vector<int> i) const
       { vector<double> yy(NYSets());
         for(int j=1;j<=static_cast<int>(NYSets());j++)
            { if(i[j-1]<1 || i[j-1]>static_cast<int>(NY(j))){ throw OUT_OF_RANGE<int>(i[j-1],1,NY(j),string("Y"),CLASS_NAME);}
              yy[j-1]=_pY(j,i[j]);
             } 
         return yy;
        }         

void YDATA_MULTIPLESETS::SetY(int i,int j,double Y) const
       { if(i<1 || i>static_cast<int>(NYSets())){ throw OUT_OF_RANGE<int>(i,1,NYSets(),string("SetY"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NY(i))){ throw OUT_OF_RANGE<int>(j,1,NY(i),string("SetY"),CLASS_NAME);}
         return _pSetY(i,j,Y);
        }

void YDATA_MULTIPLESETS::SetY(int i,const vector<double> &Y) const
       { if(i<1 || i>static_cast<int>(NYSets())){ throw OUT_OF_RANGE<int>(i,1,NYSets(),string("SetY"),CLASS_NAME);}
         return _pSetY(i,Y);
        }

// *******************************************************

void YDATA_MULTIPLESETS::AddY(int i,const vector<double> &Y)
     { vector<double> newy(_pY(i)); 
       for(int j=1;j<=static_cast<int>(Y.size());j++){ newy.push_back(Y[j-1]);}
       _pSetY(i,newy);
      }

void YDATA_MULTIPLESETS::AddY(const vector<vector<double> > &Y)
     { if(Y.size()!=NYSets()){ throw DIFFERENT_LENGTHS("AddY",CLASS_NAME,152);}
       for(int i=1;i<=static_cast<int>(NYSets());i++)
             { vector<double> newy(_pY(i));            
                for(int j=1;j<=static_cast<int>(Y[i].size());j++){ newy.push_back(Y[i-1][j-1]);}
               _pSetY(i,newy);
             }
      }

// *******************************************************

// shift, rescale the data      
void YDATA_MULTIPLESETS::ShiftY(int i,const double &S)
     { if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME,158);}
       if(Empty(i)==true){ throw EMPTY("ShiftY",CLASS_NAME,159);}

       vector<double> YT(NY(i));
       for(int j=1;j<=static_cast<int>(NY(i));j++){ YT[j-1]=_pY(i,j)+S;}
       _pSetY(i,YT);
      }

void YDATA_MULTIPLESETS::ShiftY(const double &S)
     { if(Empty()==true){ throw EMPTY("ShiftY",CLASS_NAME,167);}
       for(int i=1;i<=static_cast<int>(NYSets());i++)
          { if(Empty(i)==true){ throw EMPTY("ShiftY",CLASS_NAME,169);}
            vector<double> YT(NY(i));
            for(int j=1;j<=static_cast<int>(NY(i));j++){ YT[j-1]=_pY(i,j)+S;}
            _pSetY(i,YT);
          }
      }

void YDATA_MULTIPLESETS::RescaleY(int i,const double &S)
     { if(Empty()==true){ throw EMPTY("RescaleY",CLASS_NAME,177);}
       if(Empty(i)==true){ throw EMPTY("RescaleY",CLASS_NAME,178);}

       vector<double> YT(NY(i));
       for(int j=1;j<=static_cast<int>(NY(i));j++){ YT[j-1]=_pY(i,j)*S;}
       _pSetY(i,YT);
      }

void YDATA_MULTIPLESETS::RescaleY(const double &S)
     { if(Empty()==true){ throw EMPTY("RescaleY",CLASS_NAME,186);}
       for(int i=1;i<=static_cast<int>(NYSets());i++)
          { if(Empty(i)==true){ throw EMPTY("RescaleY",CLASS_NAME,188);}
            vector<double> YT(NY(i));
            for(int j=1;j<=static_cast<int>(NY(i));j++){ YT[j-1]=_pY(i,j)*S;}
            _pSetY(i,YT);
          }
      }

void YDATA_MULTIPLESETS::AbsY(int i)
          { for(int j=1;j<=static_cast<int>(NY(i));j++){ _pSetY(i,j,Sign(_pY(i,j))*_pY(i,j));} 
           }

void YDATA_MULTIPLESETS::AbsY(void)
          { for(int i=1;i<=static_cast<int>(NYSets());i++)
                   { for(int j=1;j<=static_cast<int>(NY(i));j++){ _pSetY(i,j,Sign(_pY(i,j))*_pY(i,j));} 
                    } 
            }

// *******************************************************
// *******************************************************
// *******************************************************

string DELTAX_SINGLESET::CLASS_NAME("DELTAX_SINGLESET");

void DELTAX_SINGLESET::XDifference(void) const
     { CreateDeltaX(static_cast<int>(NX())-1);
       for(int a=1;a<=static_cast<int>(NX())-1;a++){ _pSetDeltaX(a,_pX(a+1)-_pX(a));}
       SetXDifferenced(true);
      }

double DELTAX_SINGLESET::DeltaX(int i) const
       { if(i<1 || i>static_cast<int>(NDeltaX())){ throw OUT_OF_RANGE<int>(i,1,NDeltaX(),string("DeltaX"),CLASS_NAME);}
         return _pDeltaX(i);
        }

// *******************************************************
// *******************************************************
// *******************************************************

string DELTAX_MULTIPLESETS::CLASS_NAME("DELTAX_MULTIPLESETS");

void DELTAX_MULTIPLESETS::CreateDeltaX(size_t NDeltaXSets,std::vector<size_t> NDeltaX) const
          { SetXDifferenced(true);
            CreateDeltaX(NDeltaXSets);
            for(int i=1;i<=static_cast<int>(NDeltaXSets);i++){ deltax[i-1]=std::vector<double>(NDeltaX[i-1]);}             
           }

void DELTAX_MULTIPLESETS::DestroyDeltaX(void) const
     { if(XDifferenced()==false){ return;} 
       for(int i=1;i<=static_cast<int>(NXSets());i++){ deltax[i-1].clear();} 
       SetXDifferenced(false);
      }

void DELTAX_MULTIPLESETS::XDifference(void) const
     { CreateDeltaX(NXSets());
       for(int a=1;a<=static_cast<int>(NXSets());a++)
          { CreateDeltaXSet(a,NX(a)-1);
            for(int b=1;b<=static_cast<int>(NX(a))-1;b++){ _pSetDeltaX(a,b,_pX(a,b+1)-_pX(a,b));} 
           }
       SetXDifferenced(true);
      }

vector<size_t> DELTAX_MULTIPLESETS::NDeltaX(void) const
      { vector<size_t> n(NDeltaXSets());
        for(int i=1;i<=static_cast<int>(NDeltaXSets());i++){ n[i-1]=deltax[i-1].size();}
        return n;
       }

double DELTAX_MULTIPLESETS::DeltaX(int i,int j) const
       { if(i<1 || i>static_cast<int>(NDeltaXSets())){ throw OUT_OF_RANGE<int>(i,1,NDeltaXSets(),string("DeltaX"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NDeltaX(i))){ throw OUT_OF_RANGE<int>(j,1,NDeltaX(i),string("DeltaX"),CLASS_NAME);}
         return _pDeltaX(i,j);
        }

vector<double> DELTAX_MULTIPLESETS::DeltaX(int i) const
       { if(i<1 || i>static_cast<int>(NDeltaXSets())){ throw OUT_OF_RANGE<int>(i,1,NDeltaXSets(),string("DeltaX"),CLASS_NAME);}
         return _pDeltaX(i);
        }

vector<double> DELTAX_MULTIPLESETS::DeltaX(vector<int> i) const
       { vector<double> dx(NDeltaXSets());
         for(int j=1;j<=static_cast<int>(NDeltaXSets());j++)
            { if(i[j-1]<1 || i[j-1]>static_cast<int>(NDeltaX(j))){ throw OUT_OF_RANGE<int>(i[j-1],1,NDeltaX(j),string("DeltaX"),CLASS_NAME);}
              dx[j-1]=_pDeltaX(j,i[j]);
             } 
         return dx;
        }  

// *******************************************************
// *******************************************************
// *******************************************************

string DELTAY_SINGLESET::CLASS_NAME("DELTAY_SINGLESET");

void DELTAY_SINGLESET::YDifference(void) const
     { CreateDeltaY(static_cast<int>(NY())-1);
       for(int a=1;a<=static_cast<int>(NY())-1;a++){ _pSetDeltaY(a,_pY(a+1)-_pY(a));}
       SetYDifferenced(true);
      }

double DELTAY_SINGLESET::DeltaY(int i) const
       { if(i<1 || i>static_cast<int>(NDeltaY())){ throw OUT_OF_RANGE<int>(i,1,NDeltaY(),string("DeltaY"),CLASS_NAME);}
         return _pDeltaY(i);
        }

// *******************************************************
// *******************************************************
// *******************************************************

string DELTAY_MULTIPLESETS::CLASS_NAME("DELTAY_MULTIPLESETS");

void DELTAY_MULTIPLESETS::CreateDeltaY(size_t NDeltaYSets,std::vector<size_t> NDeltaY) const
          { CreateDeltaY(NDeltaYSets);
            for(int i=1;i<=static_cast<int>(NDeltaYSets);i++){ deltay[i-1]=std::vector<double>(NDeltaY[i-1]);} 
            SetYDifferenced(false);
           }

void DELTAY_MULTIPLESETS::DestroyDeltaY(void) const
     { if(YDifferenced()==false){ return;} 
       for(int i=1;i<=static_cast<int>(NYSets());i++){ deltay[i-1].clear();} 
       SetYDifferenced(false);
      }

void DELTAY_MULTIPLESETS::YDifference(void) const
     { CreateDeltaY(NYSets());
       for(int a=1;a<=static_cast<int>(NYSets());a++)
          { CreateDeltaYSet(a,static_cast<int>(NY(a))-1);
            for(int b=1;b<=static_cast<int>(NY(a))-1;b++){ _pSetDeltaY(a,b,_pY(a,b+1)-_pY(a,b));} 
           }
       SetYDifferenced(true);
      }

vector<size_t> DELTAY_MULTIPLESETS::NDeltaY(void) const
      { vector<size_t> n(NDeltaYSets());
        for(int i=1;i<=static_cast<int>(NDeltaYSets());i++){ n[i-1]=deltay[i-1].size();}
        return n;
       }

double DELTAY_MULTIPLESETS::DeltaY(int i,int j) const
       { if(i<1 || i>static_cast<int>(NDeltaYSets())){ throw OUT_OF_RANGE<int>(i,1,NDeltaYSets(),string("DeltaY"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NDeltaY(i))){ throw OUT_OF_RANGE<int>(j,1,NDeltaY(i),string("DeltaY"),CLASS_NAME);}
         return _pDeltaY(i,j);
        }

vector<double> DELTAY_MULTIPLESETS::DeltaY(int i) const
       { if(i<1 || i>static_cast<int>(NDeltaYSets())){ throw OUT_OF_RANGE<int>(i,1,NYSets(),string("DeltaY"),CLASS_NAME);}
         return _pDeltaY(i);
        }

vector<double> DELTAY_MULTIPLESETS::DeltaY(vector<int> i) const
       { vector<double> dy(NDeltaYSets());
         for(int j=1;j<=static_cast<int>(NDeltaYSets());j++)
            { if(i[j-1]<1 || i[j-1]>static_cast<int>(NDeltaY(j))){ throw OUT_OF_RANGE<int>(i[j-1],1,NDeltaY(j),string("DeltaY"),CLASS_NAME);}
              dy[j-1]=_pDeltaY(j,i[j]);
             } 
         return dy;
        }  

// *******************************************************
// *******************************************************
// *******************************************************

string XLIMITS_SINGLESET::CLASS_NAME("XLIMITS_SINGLESET");

void XLIMITS_SINGLESET::XLimit(void) const
       { if(Empty()==true){ throw EMPTY("XLimit");}
         double xmin=_pX(1);
         double xmax=_pX(1);
         for(int i=2;i<=static_cast<int>(NX());i++)
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

string XLIMITS_MULTIPLESETS::CLASS_NAME("XLIMITS_MULTIPLESETS");

void XLIMITS_MULTIPLESETS::XLimit(void) const
       { if(Empty()==true){ throw EMPTY("XLimit");}
         vector<double> xmin(NXSets());
         vector<double> xmax(NXSets());
         for(int i=1;i<=static_cast<int>(NXSets());i++)
                { xmin[i-1]=_pX(i,1);
                  xmax[i-1]=_pX(i,1);
                  for(int j=2;j<=static_cast<int>(NX(i));j++)
                     { if(_pX(i,j)<xmin[i-1]){ xmin[i-1]=_pX(i,j);}
                       if(_pX(i,j)>xmax[i-1]){ xmax[i-1]=_pX(i,j);}
                      }
                 }
         SetXMin(xmin);
         SetXMax(xmax);
         SetXLimited(true);
        }

double XLIMITS_MULTIPLESETS::XMin(int i) const
       { if(i<1 || i>static_cast<int>(NXSets())){ throw OUT_OF_RANGE<int>(i,1,NXSets(),string("XMin"),CLASS_NAME);}
         if(XLimited()==false){ XLimit();} 
         return _pXMin(i);
        }

double XLIMITS_MULTIPLESETS::XMax(int i) const
       { if(i<1 || i>static_cast<int>(NXSets())){ throw OUT_OF_RANGE<int>(i,1,NXSets(),string("XMin"),CLASS_NAME);}
         if(XLimited()==false){ XLimit();} 
         return _pXMax(i);
        }

vector<double> XLIMITS_MULTIPLESETS::XMin(void) const
       { if(XLimited()==false){ XLimit();} return _pXMin();}

vector<double> XLIMITS_MULTIPLESETS::XMax(void) const
       { if(XLimited()==false){ XLimit();} return _pXMax();}

// *******************************************************
// *******************************************************
// *******************************************************

string YLIMITS_SINGLESET::CLASS_NAME("YLIMITS_SINGLESET");

void YLIMITS_SINGLESET::YLimit(void) const
       { if(Empty()==true){ throw EMPTY("YLimit");}
         double ymin=_pY(1);
         double ymax=_pY(1);
         for(int i=2;i<=static_cast<int>(NY());i++)
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

string YLIMITS_MULTIPLESETS::CLASS_NAME("YLIMITS_MULTIPLESETS");

void YLIMITS_MULTIPLESETS::YLimit(void) const
       { if(Empty()==true){ throw EMPTY("YLimit");}
         vector<double> ymin(NYSets());
         vector<double> ymax(NYSets());
         for(int i=1;i<=static_cast<int>(NYSets());i++)
          { ymin[i-1]=_pY(i,1);
            ymax[i-1]=_pY(i,1);
            for(int j=2;j<=static_cast<int>(NY(i));j++)
               { if(_pY(i,j)<ymin[i-1]){ ymin[i-1]=_pY(i,j);}
                 if(_pY(i,j)>ymax[i-1]){ ymax[i-1]=_pY(i,j);}
                }
           }
         SetYMin(ymin);
         SetYMax(ymax);
         SetYLimited(true);
        }

double YLIMITS_MULTIPLESETS::YMin(int i) const
       { if(i<1 || i>static_cast<int>(NYSets())){ throw OUT_OF_RANGE<int>(i,1,NYSets(),string("YMin"),CLASS_NAME);}
         if(YLimited()==false){ YLimit();} 
         return _pYMin(i);
        }

double YLIMITS_MULTIPLESETS::YMax(int i) const
       { if(i<1 || i>static_cast<int>(NYSets())){ throw OUT_OF_RANGE<int>(i,1,NYSets(),string("YMax"),CLASS_NAME);}
         if(YLimited()==false){ YLimit();} 
         return _pYMax(i);
        }

vector<double> YLIMITS_MULTIPLESETS::YMin(void) const
       { if(YLimited()==false){ YLimit();} return _pYMin();}

vector<double> YLIMITS_MULTIPLESETS::YMax(void) const
       { if(YLimited()==false){ YLimit();} return _pYMax();}

// *******************************************************
// *******************************************************
// *******************************************************

string XYLIMITS_SINGLEXSET_SINGLEYSET::CLASS_NAME("XYLIMITS_SINGLEXSET_SINGLEYSET");

void XYLIMITS_SINGLEXSET_SINGLEYSET::XYLimit(void) const
       { double xmin=X(1), xmax=X(1), yatxmin=Y(1), yatxmax=Y(1);
         double ymin=Y(1), ymax=Y(1), xatymin=X(1), xatymax=X(1);
         for(int i=2;i<=static_cast<int>(NX());i++)
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

string XYLIMITS_SINGLEXSET_MULTIPLEYSETS::CLASS_NAME("XYLIMITS_SINGLEXSET_MULTIPLEYSETS");

void XYLIMITS_SINGLEXSET_MULTIPLEYSETS::XYLimit(void) const
       { double xmin=X(1), xmax=X(1);
         vector<double> ymin(NYSets()), ymax(NYSets()); 
         vector<double> yatxmin(NYSets()), yatxmax(NYSets());

         for(int j=1;j<=static_cast<int>(NYSets());j++)
            { ymin[j-1]=Y(j,1); 
              ymax[j-1]=Y(j,1);
              yatxmin[j-1]=Y(j,1);
              yatxmax[j-1]=Y(j,1);
             } 

         for(int i=2;i<=static_cast<int>(NX());i++)
            { if(X(i)<xmin){ xmin=X(i);  for(int j=1;j<=static_cast<int>(NYSets());j++){ yatxmin[j-1]=Y(j,i);} }
              if(X(i)>xmax){ xmax=X(i);  for(int j=1;j<=static_cast<int>(NYSets());j++){ yatxmax[j-1]=Y(j,i);} }
             }
         for(int j=2;j<=static_cast<int>(NYSets());j++)
            { for(int i=1;i<=static_cast<int>(NY(j));i++)
                 { if(Y(j,i)<ymin[j-1]){ ymin[j-1]=Y(j,i);} 
                   if(Y(j,i)>ymax[j-1]){ ymax[j-1]=Y(j,i);}
                  }
             }

         SetYAtXMin(yatxmin);
         SetYAtXMax(yatxmax);

         SetXYLimited(true);
        }

double XYLIMITS_SINGLEXSET_MULTIPLEYSETS::YAtXMin(int i) const { if(XYLimited()==false){ XYLimit();} return _pYAtXMin(i);}
vector<double> XYLIMITS_SINGLEXSET_MULTIPLEYSETS::YAtXMin(void) const { if(XYLimited()==false){ XYLimit();} return _pYAtXMin();}

double XYLIMITS_SINGLEXSET_MULTIPLEYSETS::YAtXMax(int i) const { if(XYLimited()==false){ XYLimit();} return _pYAtXMax(i);}
vector<double> XYLIMITS_SINGLEXSET_MULTIPLEYSETS::YAtXMax(void) const { if(XYLimited()==false){ XYLimit();} return _pYAtXMax();}

// *******************************************************
// *******************************************************
// *******************************************************

string GRADEBASE_SINGLEXSET_SINGLEYSET::CLASS_NAME("GRADEBASE_SINGLEXSET_SINGLEYSET");

// **************

double GRADEBASE_SINGLEXSET_SINGLEYSET::GL(int i) const
       { if(i<1 || i>static_cast<int>(NGL())){ throw OUT_OF_RANGE<int>(i,1,NGL(),string("GL"),CLASS_NAME);}
         if(Graded()==false){ Grade();}
         return _pGL(i);
        }

double GRADEBASE_SINGLEXSET_SINGLEYSET::GR(int i) const
       { if(i<1 || i>static_cast<int>(NGR())){ throw OUT_OF_RANGE<int>(i,1,NGR(),string("GR"),CLASS_NAME);}
         if(Graded()==false){ Grade();}
         return _pGR(i);
        }

double GRADEBASE_SINGLEXSET_SINGLEYSET::G(int i) const
       { if(i<1 || i>static_cast<int>(NG())){ throw OUT_OF_RANGE<int>(i,1,NG(),string("G"),CLASS_NAME);}
         if(Graded()==false){ Grade();}
         if(i==1){ return _pGL(i);}
         if(i==static_cast<int>(NG())){ return _pGR(i-1);}
         return 0.5*(_pGL(i)+_pGR(i-1));
        }

vector<double> GRADEBASE_SINGLEXSET_SINGLEYSET::GL(void) const
       { if(Graded()==false){ Grade();}
         return _pGL();
        }

vector<double> GRADEBASE_SINGLEXSET_SINGLEYSET::GR(void) const
       { if(Graded()==false){ Grade();}
         return _pGR();
        }

// **************

void GRADEBASE_SINGLEXSET_SINGLEYSET::SetGL(int i,double GL) const
       { if(i<1 || i>static_cast<int>(NGL())){ throw OUT_OF_RANGE<int>(i,1,NGL(),string("SetGL"),CLASS_NAME);}
         _pSetGL(i,GL);
        }

void GRADEBASE_SINGLEXSET_SINGLEYSET::SetGL(const vector<double> &GL) const
       { if(GL.size()!=NGL()){ throw DIFFERENT_LENGTHS("SetGL",CLASS_NAME);}
         _pSetGL(GL);
        }

void GRADEBASE_SINGLEXSET_SINGLEYSET::SetGR(int i,double GR) const
       { if(i<1 || i>static_cast<int>(NGR())){ throw OUT_OF_RANGE<int>(i,1,NGR(),string("SetGR"),CLASS_NAME);}
         _pSetGR(i,GR);
        }

void GRADEBASE_SINGLEXSET_SINGLEYSET::SetGR(const vector<double> &GR) const
       { if(GR.size()!=NGR()){ throw DIFFERENT_LENGTHS("SetGR",CLASS_NAME);}
         _pSetGR(GR);
        }


void GRADEBASE_SINGLEXSET_SINGLEYSET::SetG(int i,double G) const 
       { if(i<1 || i>static_cast<int>(NG())){ throw OUT_OF_RANGE<int>(i,1,NG(),string("SetG"),CLASS_NAME);}
         if(i==1){ _pSetGL(i,G);}
         else{ if(i==static_cast<int>(NG())){ _pSetGR(i-1,G);} 
               else{ _pSetGL(i,G); _pSetGR(i-1,G);}
              }
        }

void GRADEBASE_SINGLEXSET_SINGLEYSET::LineGrade(void) const
     { _pSetGL(1,_pDeltaY(1)/_pDeltaX(1));
       for(int i=2;i<=static_cast<int>(NG())-1;i++)
          { _pSetGL(i,_pDeltaY(i)/_pDeltaX(i));
            _pSetGR(i-1,_pDeltaY(i)/_pDeltaX(i));
           }
       _pSetGR(NG()-1,_pDeltaY(NG()-1)/_pDeltaX(NG()-1));

       SetGraded(true);
      }

void GRADEBASE_SINGLEXSET_SINGLEYSET::KochanekBartelsGrade(double t,double b,double c) const
     { double vL,vR;
       _pSetGL(1,_pDeltaY(1)/_pDeltaX(1) );       
       for(int i=2;i<=static_cast<int>(NG())-1;i++)
          { vL=_pDeltaY(i-1)/_pDeltaX(i-1);
            vR=_pDeltaY(i)/_pDeltaX(i);
            _pSetGL(i,( (1.-t)*(1.+b)*(1.-c)*_pDeltaX(i)*vL + (1.-t)*(1.-b)*(1.+c)*_pDeltaX(i-1)*vR ) / (_pDeltaX(i-1) + _pDeltaX(i)) );
            _pSetGR(i-1,( (1.-t)*(1.+b)*(1.+c)*_pDeltaX(i)*vL + (1.-t)*(1.-b)*(1.-c)*_pDeltaX(i-1)*vR ) / (_pDeltaX(i-1) + _pDeltaX(i)) );
           }
       _pSetGR(NG()-1,_pDeltaY(NG()-1)/_pDeltaX(NG()-1));

       SetGraded(true);
      }

// *******************************************************************

/*void LINEGRADE_SINGLEXSET_SINGLEYSET::Grade(void) const
     { CreateG(NX());
       int i;
       for(i=1;i<=static_cast<int>(NG())-1;i++){ _pSetG( i,_pDeltaY(i)/_pDeltaX(i) );}
       _pSetG( NG(),_pG(NG()-1) );
       
       SetGraded(true);       
      }*/

// *******************************************************

/*void LOCALGRADE_SINGLEXSET_SINGLEYSET::Grade(void) const
     { CreateG(NX());
       if(NG()<=5){ throw TOO_FEW_POINTS("Grade","LOCALGRADE_SINGLEXSET_SINGLEYSET");}

       vector<double> XX(5), YY(5), D;

       int i,xoffset;
       for(i=1;i<=static_cast<int>(NG());i++)
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
      }*/

// *******************************************************
// *******************************************************
// *******************************************************

string GRADEBASE_SINGLEXSET_MULTIPLEYSETS::CLASS_NAME("GRADEBASE_SINGLEXSET_MULTIPLEYSETS");

// ***************

void GRADEBASE_SINGLEXSET_MULTIPLEYSETS::DestroyG(void) const
     { if(Graded()==false){ return;}
       for(int i=1;i<=static_cast<int>(NGSets());i++){ gL[i-1].clear(); gR[i-1].clear();} 
       SetGraded(false);
      }

// ***************

double GRADEBASE_SINGLEXSET_MULTIPLEYSETS::GL(int i,int j) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("G"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NGL(i))){ throw OUT_OF_RANGE<int>(j,1,NGL(i),string("GL"),CLASS_NAME);}
         if(Graded()==false){ Grade();}
         return _pGL(i,j);
        }

double GRADEBASE_SINGLEXSET_MULTIPLEYSETS::GR(int i,int j) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("G"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NGR(i))){ throw OUT_OF_RANGE<int>(j,1,NGR(i),string("GR"),CLASS_NAME);}
         if(Graded()==false){ Grade();}
         return _pGR(i,j);
        }

double GRADEBASE_SINGLEXSET_MULTIPLEYSETS::G(int i,int j) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("G"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NG(i))){ throw OUT_OF_RANGE<int>(j,1,NG(i),string("G"),CLASS_NAME);}
         if(Graded()==false){ Grade();}

         if(j==1){ return _pGL(i,j);}
         if(j==static_cast<int>(NG(i))){ return _pGR(i,j-1);}
         return 0.5*(_pGL(i,j)+_pGR(i,j-1));
        }

vector<double> GRADEBASE_SINGLEXSET_MULTIPLEYSETS::GL(int i) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("G"),CLASS_NAME);}
         if(Graded()==false){ Grade();}
         return _pGL(i);
        }

vector<double> GRADEBASE_SINGLEXSET_MULTIPLEYSETS::GR(int i) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("G"),CLASS_NAME);}
         if(Graded()==false){ Grade();}
         return _pGR(i);
        }

vector<vector<double> > GRADEBASE_SINGLEXSET_MULTIPLEYSETS::GL(void) const
       { if(Graded()==false){ Grade();}
         return _pGL();
        }

vector<vector<double> > GRADEBASE_SINGLEXSET_MULTIPLEYSETS::GR(void) const
       { if(Graded()==false){ Grade();}
         return _pGR();
        }

// **************

void GRADEBASE_SINGLEXSET_MULTIPLEYSETS::SetGL(int i,int j,double GL) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("GL"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NGL(i))){ throw OUT_OF_RANGE<int>(j,1,NGL(i),string("SetGL"),CLASS_NAME);}
         _pSetGL(i,j,GL);
        }

void GRADEBASE_SINGLEXSET_MULTIPLEYSETS::SetGL(int i,const vector<double> &GL) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("GL"),CLASS_NAME);}
         if(GL.size()!=NGL(i)){ throw DIFFERENT_LENGTHS("SetGL",CLASS_NAME);}
         _pSetGL(i,GL);
        }

void GRADEBASE_SINGLEXSET_MULTIPLEYSETS::SetGR(int i,int j,double GR) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("GR"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NGR(i))){ throw OUT_OF_RANGE<int>(j,1,NGR(i),string("SetGR"),CLASS_NAME);}
         _pSetGR(i,j,GR);
        }

void GRADEBASE_SINGLEXSET_MULTIPLEYSETS::SetGR(int i,const vector<double> &GR) const
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("GR"),CLASS_NAME);}
         if(GR.size()!=NGR(i)){ throw DIFFERENT_LENGTHS("SetGR",CLASS_NAME);}
         _pSetGR(i,GR);
        }

void GRADEBASE_SINGLEXSET_MULTIPLEYSETS::SetG(int i,int j,double G) const 
       { if(i<1 || i>static_cast<int>(NGSets())){ throw OUT_OF_RANGE<int>(i,1,NGSets(),string("G"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NG(i))){ throw OUT_OF_RANGE<int>(j,1,NG(i),string("SetG"),CLASS_NAME);}
         if(j==1){ _pSetGL(i,j,G);}
         else{ if(j==static_cast<int>(NG(i))){ _pSetGR(i,j-1,G);} 
               else{ _pSetGL(i,j,G); _pSetGR(i,j-1,G);}
              }
        }

void GRADEBASE_SINGLEXSET_MULTIPLEYSETS::LineGrade(void) const
     { for(int j=1;j<=static_cast<int>(NGSets());j++){ 
           _pSetGL(j,1,_pDeltaY(j,1)/_pDeltaX(1));
           for(int i=2;i<=static_cast<int>(NG(j))-1;i++)
              { _pSetGL(j,i,_pDeltaY(j,i)/_pDeltaX(i));
                _pSetGR(j,i-1,_pDeltaY(j,i)/_pDeltaX(i));
               }
           _pSetGR(j,NG(j)-1,_pDeltaY(j,NG(j)-1)/_pDeltaX(NG(j)-1));
          }

       SetGraded(true);
      }

void GRADEBASE_SINGLEXSET_MULTIPLEYSETS::KochanekBartelsGrade(double t,double b,double c) const
     { double vL,vR;
       for(int j=1;j<=static_cast<int>(NGSets());j++){ 
           _pSetGL(j,1,_pDeltaY(j,1)/_pDeltaX(1) );       
           for(int i=2;i<=static_cast<int>(NG(j))-1;i++)
              { vL=_pDeltaY(j,i-1)/_pDeltaX(i-1);
                vR=_pDeltaY(j,i)/_pDeltaX(i);
                _pSetGL(j,i,( (1.-t)*(1.+b)*(1.-c)*_pDeltaX(i)*vL + (1.-t)*(1.-b)*(1.+c)*_pDeltaX(i-1)*vR ) / (_pDeltaX(i-1) + _pDeltaX(i)) );
                _pSetGR(j,i-1,( (1.-t)*(1.+b)*(1.+c)*_pDeltaX(i)*vL + (1.-t)*(1.-b)*(1.-c)*_pDeltaX(i-1)*vR ) / (_pDeltaX(i-1) + _pDeltaX(i)) );
               }
           _pSetGR(j,NG(j)-1,_pDeltaY(j,NG(j)-1)/_pDeltaX(NG(j)-1));
          } 

       SetGraded(true);
      }

// *******************************************************

/*void LINEGRADE_SINGLEXSET_MULTIPLEYSETS::Grade(void) const
     { CreateG(NYSets());
       for(int i=1;i<=NGSets();i++)
          { CreateGSet(i,NY(i));
             int j;
                    for(j=1;j<=NG(i)-1;j++){ _pSetG( i,j,_pDeltaY(i,j)/_pDeltaX(j) );}
             _pSetG( i,NG(i),_pG(i,NG(i)-1) ); 
           }
       
       SetGraded(true);       
      }*/

// *******************************************************

/*void LOCALGRADE_SINGLEXSET_MULTIPLEYSETS::Grade(void) const
     { CreateG(NYSets());

       vector<double> XX(5), YY(5), D;

       for(int i=1;i<=NGSets();i++)
          { CreateGSet(i,NY(i));

             if(NG(i)<=5){ throw TOO_FEW_POINTS("Grade","LOCALGRADE_SINGLEXSET_MULTIPLEYSETS");}

             int j,xoffset;
                    for(j=1;j<=NG(i);j++)
                       { if(j==1){ xoffset=0;} 
                         else{ if(j==2){ xoffset=-1;}  
                               else{ if(j==NG(i)-1){ xoffset=-3;} 
                                     else{ if(j==NG(i)){ xoffset=-4;} 
                                           else{ xoffset=-2;} 
                                          } 
                                    } 
                              }
                         
                          for(int a=0;a<=4;a++){ XX[a]=_pX(j+a+xoffset)-_pX(j); YY[a]=_pY(i,j+a+xoffset)-_pY(i,j);}
                          D=FiniteDifference1D(0.,XX,YY);
                          _pSetG(i,j,D[1]);
                         }
            }

       SetGraded(true);
      }*/

// *******************************************************
// *******************************************************
// *******************************************************

string SPLINE_SINGLEXSET_SINGLEYSET::CLASS_NAME("SPLINE_SINGLEXSET_SINGLEYSET");

void SPLINE_SINGLEXSET_SINGLEYSET::DestroyA(void) const
      { if(Fitted()==false){ return;}
        for(int i=1;i<=static_cast<int>(NA());i++){ a[i-1].clear();} 
        SetFitted(false);
       }

// *************

double SPLINE_SINGLEXSET_SINGLEYSET::Interpolate(int i,double T) const
       { double I=A(i,NParameters(i)-1);
         for(int j=static_cast<int>(NParameters(i))-2;j>=0;j--){ (I*=T)+=A(i,j);}
         return I;
        }

double SPLINE_SINGLEXSET_SINGLEYSET::Derivative(int i,double T) const
       { double D=A(i,NParameters(i)-1)*(NParameters(i)-1.);
         for(int j=static_cast<int>(NParameters(i))-2;j>=1;j--){ (D*=T)+=A(i,j)*j;}
         return D;
        }

double SPLINE_SINGLEXSET_SINGLEYSET::SecondDerivative(int i,double T) const
       { double SD=A(i,NParameters(i)-1)*(NParameters(i)-1.)*(NParameters(i)-2.);
         for(int j=static_cast<int>(NParameters(i))-2;j>=2;j--){ (SD*=T)+=A(i,j)*j*(j-1.);}
         return SD;
        }

double SPLINE_SINGLEXSET_SINGLEYSET::Integral(int i,double T0,double T1) const
       { double I1,I0;
         I0=I1=A(i,NParameters(i)-1)/(NParameters(i)-1);
         for(int j=static_cast<int>(NParameters(i))-2;j>=0;j--){ (I0*=T0)+=A(i,j)/(j+1.); (I1*=T1)+=A(i,j)/(j+1.);}
         I0*=T0; I1*=T1;
         return I1-I0;
        }

double SPLINE_SINGLEXSET_SINGLEYSET::A(int i,int j) const 
       { if(i<1 || i>static_cast<int>(NA())){ throw OUT_OF_RANGE<int>(i,1,NA(),string("A"),CLASS_NAME);}
         if(j<0 || j>static_cast<int>(NParameters(i))-1){ throw OUT_OF_RANGE<int>(j,0,NParameters(i)-1,string("A"),CLASS_NAME);}
         if(Fitted()==false){ Fit();} 
         return _pA(i,j);
        }

vector<double> SPLINE_SINGLEXSET_SINGLEYSET::A(int i) const
       { if(i<1 || i>static_cast<int>(NA())){ throw OUT_OF_RANGE<int>(i,1,NA(),string("A"),CLASS_NAME);}
         if(Fitted()==false){ Fit();} 
         return _pA(i);
        }

vector<vector<double> > SPLINE_SINGLEXSET_SINGLEYSET::A(void) const
       { if(Fitted()==false){ Fit();} 
         return _pA();
        }

// *******************************************************
// *******************************************************
// *******************************************************

string SPLINE_SINGLEXSET_MULTIPLEYSETS::CLASS_NAME("SPLINE_SINGLEXSET_MULTIPLEYSETS");

void SPLINE_SINGLEXSET_MULTIPLEYSETS::DestroyA(void) const
     { if(Fitted()==false){ return ;}
       for(int i=1;i<=static_cast<int>(NASets());i++)
             { for(int j=1;j<=static_cast<int>(NA(i));j++){ a[i-1][j-1].clear();} 
              } 
       SetFitted(false);
      }

// *************

double SPLINE_SINGLEXSET_MULTIPLEYSETS::Interpolate(int j,int i,double T) const
       { double I=A(j,i,NParameters(j,i)-1);
         for(int k=static_cast<int>(NParameters(j,i))-2;k>=0;k--){ (I*=T)+=A(j,i,k);}
         return I;
        }

vector<double> SPLINE_SINGLEXSET_MULTIPLEYSETS::Interpolate(int i,double T) const
       { vector<double> I(NASets());
          int j;
          for(j=1;j<=static_cast<int>(NASets());j++)
                { I[j-1]=A(j,i,NParameters(j,i)-1);
                   for(int k=static_cast<int>(NParameters(j,i))-2;k>=0;k--){ (I[j-1]*=T)+=A(j,i,k);}
                 }
         return I;
        }

double SPLINE_SINGLEXSET_MULTIPLEYSETS::Derivative(int j,int i,double T) const
       { double D=A(j,i,NParameters(j,i)-1)*(NParameters(j,i)-1.);
          for(int k=static_cast<int>(NParameters(j,i))-2;k>=1;k--){ (D*=T)+=A(j,i,k)*k;}
          return D;
         }

vector<double> SPLINE_SINGLEXSET_MULTIPLEYSETS::Derivative(int i,double T) const
       { vector<double> D(NASets());
          for(int j=1;j<=static_cast<int>(NASets());j++)
                { D[j-1]=A(j,i,NParameters(j,i)-1)*(NParameters(j,i)-1.);
                   for(int k=static_cast<int>(NParameters(j,i))-2;k>=1;k--){ (D[j-1]*=T)+=A(j,i,k)*k;}
                  }
         return D;
        }

double SPLINE_SINGLEXSET_MULTIPLEYSETS::SecondDerivative(int j,int i,double T) const
       { double SD=A(j,i,NParameters(j,i)-1)*(NParameters(j,i)-1.)*(NParameters(j,i)-2.);
          for(int k=static_cast<int>(NParameters(j,i))-2;k>=2;k--){ (SD*=T)+=A(j,i,k)*k*(k-1.);}
          return SD;
         }

vector<double> SPLINE_SINGLEXSET_MULTIPLEYSETS::SecondDerivative(int i,double T) const
       { vector<double> SD(NASets());
          for(int j=1;j<=static_cast<int>(NASets());j++)
                { SD[j-1]=A(j,i,NParameters(j,i)-1)*(NParameters(j,i)-1.)*(NParameters(j,i)-2.);
                   for(int k=static_cast<int>(NParameters(j,i))-2;k>=2;k--){ (SD[j-1]*=T)+=A(j,i,k)*k*(k-1.);}
                 }
         return SD;
        }

double SPLINE_SINGLEXSET_MULTIPLEYSETS::Integral(int j,int i,double T0,double T1) const
       { double I1,I0;
         I0=I1=A(j,i,NParameters(j,i)-1)/(NParameters(j,i)-1.);
         for(int k=static_cast<int>(NParameters(j,i))-2;k>=0;k--){ (I0*=T0)+=A(j,i,k)/(k+1.); (I1*=T1)+=A(j,i,k)/(k+1.);}
         I0*=T0; I1*=T1;
         return I1-I0;
        }

vector<double> SPLINE_SINGLEXSET_MULTIPLEYSETS::Integral(int i,double T0,double T1) const
       { double I1,I0;
         vector<double> I(NASets());
         for(int j=1;j<=static_cast<int>(NASets());j++)
            { I0=I1=A(j,i,NParameters(j,i)-1)/(NParameters(j,i)-1.);
              for(int k=static_cast<int>(NParameters(j,i))-2;k>=0;k--){ (I0*=T0)+=A(j,i,k)/(k+1.); (I1*=T1)+=A(j,i,k)/(k+1.);}
              I0*=T0; I1*=T1;
              I[j-1]=I1-I0;
             }
         return I;
        }

double SPLINE_SINGLEXSET_MULTIPLEYSETS::A(int i,int j,int k) const
       { if(i<1 || i>static_cast<int>(NASets())){ throw OUT_OF_RANGE<int>(i,1,NASets(),string("A"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NA(i))){ throw OUT_OF_RANGE<int>(j,1,NA(i),string("A"),CLASS_NAME);}
         if(k<0 || k>static_cast<int>(NParameters(i,j))-1){ throw OUT_OF_RANGE<int>(k,0,NParameters(i,j)-1,string("A"),CLASS_NAME);}         
         if(Fitted()==false){ Fit();} 
         return _pA(i,j,k);
        }

vector<double> SPLINE_SINGLEXSET_MULTIPLEYSETS::A(int i,int j) const
       { if(i<1 || i>static_cast<int>(NASets())){ throw OUT_OF_RANGE<int>(i,1,NASets(),string("A"),CLASS_NAME);}
         if(j<1 || j>static_cast<int>(NA(i))){ throw OUT_OF_RANGE<int>(j,1,NA(i),string("A"),CLASS_NAME);}
         if(Fitted()==false){ Fit();} 
         return _pA(i,j);
        }


vector<vector<double> > SPLINE_SINGLEXSET_MULTIPLEYSETS::A(int i) const
       { if(i<1 || i>static_cast<int>(NASets())){ throw OUT_OF_RANGE<int>(i,1,NASets(),string("A"),CLASS_NAME);}
         if(Fitted()==false){ Fit();} 
         return _pA(i);
        }

vector<vector<vector<double> > > SPLINE_SINGLEXSET_MULTIPLEYSETS::A(void) const
       { if(Fitted()==false){ Fit();} 
         return _pA();
        }

/*vector<double> SPLINE_SINGLEXSET_MULTIPLEYSETS::FindRoots(int j,int i,double Y) const
       { vector<double> Rji,R;

         if(NParameters(j,i)==2){ R = RealLinearRoots(A(j,i,1),A(j,i,0),Y);}
         else{ R = RealCubicRoots(A(j,i,3),A(j,i,2),A(j,i,1),A(j,i,0),Y);}

         for(int k=0;k<=static_cast<int>(Rji.size())-1;k++)
            { if(i<=N()-2 && _pX(i)<=R[0] && R[0]<_pX(i+1)){ Rji.push_back(Rji[k]);} 
              if(i==N()-1 && _pX(i)<=R[0] && R[0]<=_pX(i+1)){ Rji.push_back(Rji[k]);} 
             }

         return Rji;
        }

vector<vector<double> > SPLINE_SINGLEXSET_MULTIPLEYSETS::FindRoots(int i,vector<double> Y) const
       { vector<vector<double> > Ri(NASets());

         for(int j=1;j<=static_cast<int>(NASets());j++)
            { Rji = FindRoots(j,i,Y[j-1];}
              for(int k=0;k<=static_cast<int>(Rji.size())-1;k++)
                 { Ri[j].push_back(Rji[k]);}
             }

         return Ri;   
        }*/

// *******************************************************************
// *******************************************************************
// *******************************************************************

void SORT_SINGLEXSET_SINGLEYSET::Sort(void)
     { double tmpx, tmpy;
       for(int j=2;j<=static_cast<int>(NX());j++)
          { tmpx=_pX(j); tmpy=_pY(j);
            int i=j-1;
            while(i>=1 && _pX(i)>tmpx){ _pSetX(i+1,_pX(i)); _pSetY(i+1,_pY(i)); i--;};
            _pSetX(i+1,tmpx); _pSetY(i+1,tmpy);       
           }
       SetSorted(true); 
      }

void SORT_SINGLEXSET_MULTIPLEYSETS::Sort(void)
     { double tmpx;
       vector<double> tmpy(NYSets());
       for(int j=2;j<=static_cast<int>(NX());j++)
             { tmpx=_pX(j); 
       for(int k=1;k<=static_cast<int>(NYSets());k++){ tmpy[k-1]=_pY(k,j);}
       int i=j-1;
       while(i>=1 && _pX(i)>tmpx)
            { _pSetX(i+1,_pX(i)); 
              for(int k=1;k<=static_cast<int>(NYSets());k++){ _pSetY(i+1,k,_pY(k,i));}
              i--;
             };
       _pSetX(i+1,tmpx); 
       for(int k=1;k<=static_cast<int>(NYSets());k++){ _pSetY(k,i+1,tmpy[k-1]);}
           }
       SetSorted(true);
      }

// *******************************************************************
// *******************************************************************
// *******************************************************************

void READWRITE_SINGLEXSET_SINGLEYSET::Read(string filename,size_t NI)
       { ifstream fin(filename.c_str());

         if(fin){ size_t n, offset;
                  
                  double d; 
                  vector<double> data; 
                  string line;            
             
                  LOADING("Read",filename);
                  for(int a=0;a<=static_cast<int>(NI)-1;a++){ getline(fin,line);}

                  while(fin>>d){ data.push_back(d);};
                  fin.close();
             
                  if(Odd(data.size())==true){ n=static_cast<size_t>(data[0]); offset=1;} else{ n=data.size()/2; offset=0;}
                  CreateX(n);
                  CreateY(n);
  		  for(int a=0;a<=static_cast<int>(n)-1;a++){ _pSetX(a+1,data[2*a+offset]); _pSetY(a+1,data[2*a+offset+1]);}
                 }
         else{ throw CANNOT_FIND("Read",filename);}
        }

void READWRITE_SINGLEXSET_SINGLEYSET::Read(string filename,char ignore)
       { ifstream fin(filename.c_str());

         if(fin){ size_t n, offset;
                  double d; 
                  vector<double> data;             
                  string line;
             
                  LOADING("Read",filename);    
                  while(fin.peek()==ignore){ getline(fin,line);};
                  while(fin>>d){ data.push_back(d);};
                  fin.close();
             
                  if(Odd(data.size())==true){ n=static_cast<size_t>(data[0]); offset=1;} else{ n=data.size()/2; offset=0;}
                  CreateX(n);
                  CreateY(n);
  		  for(int a=0;a<=static_cast<int>(n)-1;a++){ _pSetX(a+1,data[2*a+offset]); _pSetY(a+1,data[2*a+offset+1]);}
                 }
         else{ throw CANNOT_FIND("Read",filename);}
        }

void READWRITE_SINGLEXSET_SINGLEYSET::Write(string filename)
     { ofstream fout(filename.c_str());
       if(fout){ WRITING("Write",filename); 
                     fout<<NX()<<"\n";
                     for(int i=1;i<=static_cast<int>(NX());i++){ fout<<_pX(i)<<"\t"<<_pY(i)<<"\n";}
                     fout.flush();
                     fout.close();
                   }
       else{ throw CANNOT_WRITE("Write",filename);}
      } 

// output the x and y's in point format
ostream& READWRITE_SINGLEXSET_SINGLEYSET::operator<<(ostream &os)
         { os<<NX()<<flush; os.precision(15);
           for(int a=1;a<=static_cast<int>(NX());a++){ os<<"\n"<<_pX(a)<<"\t"<<_pY(a)<<flush;}
           return os;
          }

istream& READWRITE_SINGLEXSET_SINGLEYSET::operator>>(istream &is)
         { size_t n, offset;
           double d; 
           vector<double> data;             
             
           while(is>>d){ data.push_back(d);};
             
           if(Odd(data.size())==true){ n=static_cast<size_t>(data[0]); offset=1;} else{ n=data.size()/2; offset=0;}
           CreateX(n);
           CreateY(n);

           for(int a=0;a<=static_cast<int>(n)-1;a++){ _pSetX(a+1,data[2*a+offset]); _pSetY(a+1,data[2*a+offset+1]);}

           return is;
          }

// *******************************************************
// *******************************************************
// *******************************************************

void READWRITE_SINGLEXSET_MULTIPLEYSETS::Read(string filename,size_t NI,bool read_size)
       { ifstream fin(filename.c_str());

         if(fin){ size_t n,nysets;
                  double d; 
                  string line;            
               
                  LOADING("Read",filename);
                  for(int a=0;a<=static_cast<int>(NI)-1;a++){ getline(fin,line);}

                  if(read_size == true)
                    { fin>>n>>nysets;
                      CreateX(n);
                      CreateY(nysets,n);
                     }
                  else{ n = NX(); nysets = NYSets();} 

                  for(int a=0;a<=static_cast<int>(n)-1;a++){ 
                      fin>>d;
                      _pSetX(a+1,d);
                      for(int b=0;b<=static_cast<int>(nysets)-1;b++){ fin>>d; _pSetY(b+1,a+1,d);}
                     }
                  fin.close();
                 }
         else{ throw CANNOT_FIND("Read",filename);}
       }

void READWRITE_SINGLEXSET_MULTIPLEYSETS::Read(string filename,char ignore,bool read_size)
       { ifstream fin(filename.c_str());

         if(fin){ size_t n,nysets;
                  double d; 
                  string line;        
             
                  LOADING("Read",filename);
                  while(fin.peek()==ignore){ getline(fin,line);};

                  if(read_size == true)
                    { fin>>n>>nysets;
                      CreateX(n);
                      CreateY(nysets,n);
                     }
                  else{ n = NX(); nysets = NYSets();}

                  for(int a=0;a<=static_cast<int>(n)-1;a++){ 
                      fin>>d;
                      _pSetX(a+1,d);
                      for(int b=0;b<=static_cast<int>(nysets)-1;b++){ fin>>d; _pSetY(b+1,a+1,d);}
                     }
                  fin.close();
                  TRACE("Finished Reading");
                 }
         else{ throw CANNOT_FIND("Read",filename);}
       }

void READWRITE_SINGLEXSET_MULTIPLEYSETS::Write(string filename)
     { ofstream fout(filename.c_str());
       if(fout){ WRITING("Write",filename); 
                     fout<<NX()<<NYSets();
                     for(int i=1;i<=static_cast<int>(NX());i++)
                           { fout<<"\n"<<_pX(i);
                              for(int j=1;j<=static_cast<int>(NYSets());j++){ fout<<"\t"<<_pY(j,i);}
                            }
                     fout.flush();
                     fout.close();
                   }
       else{ throw CANNOT_WRITE("Write",filename);}
      } 

// output the x and y's in point format
ostream& READWRITE_SINGLEXSET_MULTIPLEYSETS::operator<<(ostream &os)
         { os.precision(15);
           os<<NX()<<"\t"<<NYSets()<<flush;  
           for(int a=1;a<=static_cast<int>(NX());a++)
                 { os<<"\n"<<_pX(a);
                   for(int b=1;b<=static_cast<int>(NYSets());b++){ os<<"\t"<<_pY(b,a);}
                   os<<flush;
                  }
           return os;
          }

// *******************************************************************
// *******************************************************************
// *******************************************************************

void READWRITE_MULTIPLEXSETS_SINGLEYSET<2>::Read(string filename,size_t NI,bool read_size)
       { ifstream fin(filename.c_str());

         if(fin){ vector<size_t> n(2);
                  vector<int> i(2);                  
                  vector<double> dataX(2);
                  double dataY;
                  string line;  
             
                  LOADING("Read",filename);
                  for(int a=0;a<=static_cast<int>(NI)-1;a++){ getline(fin,line);}

                  if(read_size==true)
                    { fin>>n[0]>>n[1]; 
                      CreateX(2,n);
                      CreateY(n);
                     }
                  else{ n[0]=NX(1); n[1]=NX(2);}

  		  for(i[0]=1;i[0]<=static_cast<int>(n[0]);i[0]++)
                     { for(i[1]=1;i[1]<=static_cast<int>(n[1]);i[1]++)
                          { for(int a=0;a<=1;a++){ fin>>dataX[a];}
                            fin>>dataY;
                            _pSetX(1,i[0],dataX[0]);
                            _pSetX(2,i[1],dataX[1]); 
                            _pSetY(i,dataY);
                           }
                      }
                  fin.close();
                 }
         else{ throw CANNOT_FIND("Read",filename);}
        }

void READWRITE_MULTIPLEXSETS_SINGLEYSET<2>::Read(string filename,char ignore,bool read_size)
       { ifstream fin(filename.c_str());

         if(fin){ vector<size_t> n(2);
                  vector<int> i(2);         
                  vector<double> dataX(2);             
                  double dataY;
                  string line;
             
                  LOADING("Read",filename);    
                  while(fin.peek()==ignore){ getline(fin,line);};

                  if(read_size==true)
                    { fin>>n[0]>>n[1];
                      CreateX(2,n);
                      CreateY(n);
                     }
                  else{ n[0]=NX(1); n[1]=NX(2);}

  		  for(i[0]=1;i[0]<=static_cast<int>(n[0]);i[0]++)
                     { for(i[1]=1;i[1]<=static_cast<int>(n[1]);i[1]++)
                          { for(int a=0;a<=1;a++){ fin>>dataX[a];}
                            fin>>dataY;
                            _pSetX(1,i[0],dataX[0]);
                            _pSetX(2,i[1],dataX[1]); 
                            _pSetY(i,dataY);
                           }
                      }
                  fin.close();
                 }
         else{ throw CANNOT_FIND("Read",filename);}
        }

void READWRITE_MULTIPLEXSETS_SINGLEYSET<2>::Write(string filename)
     { ofstream fout(filename.c_str());
       if(fout){ WRITING("Write",filename); 

                 vector<size_t> n=NX();
                 vector<int> i(2);     
                 vector<vector<double> > dataX=X();
                 TARRAY<double,2> dataY=Y();

                 fout<<n[0]<<"\t"<<n[1];
                 
                 for(i[0]=0;i[0]<=static_cast<int>(n[0])-1;i[0]++)
                     { for(i[1]=0;i[1]<=static_cast<int>(n[1])-1;i[1]++)
                          { for(int a=0;a<=1;a++)
                               { fout<<"\n"<<dataX[0][i[0]]<<"\t"<<dataX[1][i[1]]<<"\t"<<dataY[i[0]][i[1]];}
                           }
                      }
                 fout.flush();
                 fout.close();
                }
       else{ throw CANNOT_WRITE("Write",filename);}
      } 

// output the x and y's in point format
ostream& READWRITE_MULTIPLEXSETS_SINGLEYSET<2>::operator<<(ostream &os)
         { vector<size_t> n=NX();
           vector<int> i(2);     
           vector<vector<double> > dataX=X();
           TARRAY<double,2> dataY=Y();

           os<<n[0]<<"\t"<<n[1];
                 
           for(i[0]=0;i[0]<=static_cast<int>(n[0])-1;i[0]++)
              { for(i[1]=0;i[1]<=static_cast<int>(n[1])-1;i[1]++)
                   { for(int a=0;a<=1;a++)
                        { os<<"\n"<<dataX[0][i[0]]<<"\t"<<dataX[1][i[1]]<<"\t"<<dataY[i[0]][i[1]];}
                    }
               }

           return os;
          }

istream& READWRITE_MULTIPLEXSETS_SINGLEYSET<2>::operator>>(istream &is)
         { vector<size_t> n(2);
           vector<int> i(2);
           vector<double> dataX(2);
           double dataY; 
             
           is>>n[0]>>n[1];
           CreateX(2,n);
           CreateY(n);
  	   for(i[0]=1;i[0]<=static_cast<int>(n[0]);i[0]++)
              { for(i[1]=1;i[1]<=static_cast<int>(n[1]);i[1]++)
                   { for(int a=0;a<=1;a++){ is>>dataX[a];}
                     is>>dataY;
                     _pSetX(1,i[0],dataX[0]);
                     _pSetX(2,i[1],dataX[1]); 
                     _pSetY(i,dataY);
                    }
               }

           return is;
          }

// *******************************************************************

void READWRITE_MULTIPLEXSETS_SINGLEYSET<3>::Read(string filename,size_t NI,bool read_size)
       { ifstream fin(filename.c_str());

         if(fin){ vector<size_t> n(3);
                  vector<int> i(3);                  
                  vector<double> dataX(3);
                  double dataY;
                  string line;  
             
                  LOADING("Read",filename);
                  for(int a=0;a<=static_cast<int>(NI)-1;a++){ getline(fin,line);}

                  if(read_size==true)
                    { fin>>n[0]>>n[1]>>n[2]; 
                      CreateX(3,n);
                      CreateY(n);
                     }
                  else{ n[0]=NX(1); n[1]=NX(2); n[2]=NX(3);}

  		  for(i[0]=1;i[0]<=static_cast<int>(n[0]);i[0]++)
                     { for(i[1]=1;i[1]<=static_cast<int>(n[1]);i[1]++)
                          { for(i[2]=1;i[2]<=static_cast<int>(n[2]);i[2]++)
                               { for(int a=0;a<=2;a++){ fin>>dataX[a];}
                                 fin>>dataY;
                                 _pSetX(1,i[0],dataX[0]);
                                 _pSetX(2,i[1],dataX[1]); 
                                 _pSetX(3,i[2],dataX[2]);
                                 _pSetY(i,dataY);
                                }
                           }
                      }
                  fin.close();
                 }
         else{ throw CANNOT_FIND("Read",filename);}
        }

void READWRITE_MULTIPLEXSETS_SINGLEYSET<3>::Read(string filename,char ignore,bool read_size)
       { ifstream fin(filename.c_str());

         if(fin){ vector<size_t> n(3);
                  vector<int> i(3);         
                  vector<double> dataX(3);             
                  double dataY;
                  string line;
             
                  LOADING("Read",filename);    
                  while(fin.peek()==ignore){ getline(fin,line);};

                  if(read_size==true)
                    { fin>>n[0]>>n[1]>>n[2];
                      CreateX(3,n);
                      CreateY(n);
                     }
                  else{ n[0]=NX(1); n[1]=NX(2); n[2]=NX(3);}

  		  for(i[0]=1;i[0]<=static_cast<int>(n[0]);i[0]++)
                     { for(i[1]=1;i[1]<=static_cast<int>(n[1]);i[1]++)
                          { for(i[2]=1;i[2]<=static_cast<int>(n[2]);i[2]++)
                               { for(int a=0;a<=2;a++){ fin>>dataX[a];}
                                 fin>>dataY;
                                 _pSetX(1,i[0],dataX[0]);
                                 _pSetX(2,i[1],dataX[1]); 
                                 _pSetX(3,i[2],dataX[2]);
                                 _pSetY(i,dataY);
                                }
                           }
                      }
                  fin.close();
                 }
         else{ throw CANNOT_FIND("Read",filename);}
        }

void READWRITE_MULTIPLEXSETS_SINGLEYSET<3>::Write(string filename)
     { ofstream fout(filename.c_str());
       if(fout){ WRITING("Write",filename); 

                 vector<size_t> n=NX();
                 vector<int> i(3);     
                 vector<vector<double> > dataX=X();
                 TARRAY<double,3> dataY=Y();

                 fout<<n[0]<<"\t"<<n[1]<<"\t"<<n[2];
                 
                 for(i[0]=0;i[0]<=static_cast<int>(n[0])-1;i[0]++)
                     { for(i[1]=0;i[1]<=static_cast<int>(n[1])-1;i[1]++)
                          { for(i[2]=0;i[2]<=static_cast<int>(n[2])-1;i[2]++)
                               { for(int a=0;a<=2;a++)
                                    { fout<<"\n"<<dataX[0][i[0]]<<"\t"<<dataX[1][i[1]]<<dataX[2][i[2]]<<"\t"<<dataY[i[0]][i[1]][i[2]];}
                                }
                           }
                      }
                 fout.flush();
                 fout.close();
                }
       else{ throw CANNOT_WRITE("Write",filename);}
      } 

// output the x and y's in point format
ostream& READWRITE_MULTIPLEXSETS_SINGLEYSET<3>::operator<<(ostream &os)
         { vector<size_t> n=NX();
           vector<int> i(3);     
           vector<vector<double> > dataX=X();
           TARRAY<double,3> dataY=Y();

           os<<n[0]<<"\t"<<n[1]<<"\t"<<n[2];
                 
           for(i[0]=0;i[0]<=static_cast<int>(n[0])-1;i[0]++)
              { for(i[1]=0;i[1]<=static_cast<int>(n[1])-1;i[1]++)
                   { for(i[2]=0;i[2]<=static_cast<int>(n[2])-1;i[2]++)
                        { for(int a=0;a<=2;a++)
                             { os<<"\n"<<dataX[0][i[0]]<<"\t"<<dataX[1][i[1]]<<dataX[2][i[2]]<<"\t"<<dataY[i[0]][i[1]][i[2]];}
                         }
                    }
               }

           return os;
          }

istream& READWRITE_MULTIPLEXSETS_SINGLEYSET<3>::operator>>(istream &is)
         { vector<size_t> n(3);
           vector<int> i(3);
           vector<double> dataX(3);
           double dataY; 
             
           is>>n[0]>>n[1]>>n[2];
           CreateX(3,n);
           CreateY(n);
  	   for(i[0]=1;i[0]<=static_cast<int>(n[0]);i[0]++)
              { for(i[1]=1;i[1]<=static_cast<int>(n[1]);i[1]++)
                   { for(i[2]=1;i[2]<=static_cast<int>(n[2]);i[2]++)
                        { for(int a=0;a<=2;a++){ is>>dataX[a];}
                          is>>dataY;
                          _pSetX(1,i[0],dataX[0]);
                          _pSetX(2,i[1],dataX[1]); 
                          _pSetX(3,i[2],dataX[2]);  
                          _pSetY(i,dataY);
                         }
                    }
               }

           return is;
          }

// *******************************************************************
// *******************************************************************
// *******************************************************************

int XINTERVAL_SINGLESET::XInterval(double X) const
      { if(X<XMin() || X>XMax()){ throw OUT_OF_RANGE<double>(X,XMin(),XMax(),string("XInterval"));}

        if(interval<1 || interval>static_cast<int>(NX())-1){ interval=(static_cast<int>(NX())-1)/2;}

        if(interval!=static_cast<int>(NX())-1 && X>=_pX(interval) && _pX(interval+1)>X){ return interval;}
        if(interval==static_cast<int>(NX())-1 && X>=_pX(interval) && _pX(interval+1)>=X){ return interval;}

        if(interval>=2){ if(X>=_pX(interval-1) && _pX(interval)>X){ --interval; return interval;} }
        if(interval<=static_cast<int>(NX())-3){ if(X>=_pX(interval+1) && _pX(interval+2)>X){ ++interval; return interval;} }
        if(interval==static_cast<int>(NX())-2){ if(X>=_pX(interval+1) && _pX(interval+2)>=X){ ++interval; return interval;} }
  
        if(Equality(X,_pX(1))==true){ return interval=1;} if(Equality(X,_pX(NX()))==true){ return interval=static_cast<int>(NX())-1;}

        int lower=1, upper=static_cast<int>(NX())-1;
        interval=( upper+lower + (upper-lower)%2 )/2;

        while(lower!=upper && (X<_pX(interval) || X>=_pX(interval+1)) )
             { if(X>=_pX(interval+1)){ lower=interval+1;} else{ upper=interval-1;}
               interval=( upper+lower + (upper-lower)%2 )/2;
              }

        return interval;
       }

// *******************************************************************

std::vector<int> XINTERVAL_MULTIPLESETS::XInterval(std::vector<double> X) const
      { bool found; 
        for(int i=1;i<=static_cast<int>(NXSets());i++)
           { found=false;
             if(X[i-1]<XMin(i) || X[i-1]>XMax(i)){ throw OUT_OF_RANGE<double>(X[i-1],XMin(i),XMax(i),string("XInterval"));}

             if(interval[i-1]<1 || interval[i-1]>static_cast<int>(NX(i))-1){ interval[i-1]=(static_cast<int>(NX(i))-1)/2;}

             if(interval[i-1]>=2){ if(X[i-1]>=_pX(i,interval[i-1]-1) && _pX(i,interval[i-1])>X[i-1]){ --interval[i-1]; found=true;} }
             if(interval[i-1]<=static_cast<int>(NX(i))-3){ if(X[i-1]>=_pX(i,interval[i-1]+1) && _pX(i,interval[i-1]+2)>X[i-1]){ ++interval[i-1]; found=true;} }
             if(interval[i-1]==static_cast<int>(NX(i))-2){ if(X[i-1]>=_pX(i,interval[i-1]+1) && _pX(i,interval[i-1]+2)>=X[i-1]){ ++interval[i-1]; found=true;} }
  
             if(Equality(X[i-1],_pX(i,1))==true){ interval[i-1]=1; found=true;} 
             if(Equality(X[i-1],_pX(i,NX(i)))==true){ interval[i-1]=static_cast<int>(NX(i))-1; found=true;}

             if(found==false)
               { int lower=1, upper=static_cast<int>(NX(i))-1;
                 interval[i-1]=( upper+lower + (upper-lower)%2 )/2;

                 while(lower!=upper && (X[i-1]<_pX(i,interval[i-1]) || X[i-1]>=_pX(i,interval[i-1]+1)) )
                      { if(X[i-1]>=_pX(i,interval[i-1]+1)){ lower=interval[i-1]+1;} else{ upper=interval[i-1]-1;}
                        interval[i-1]=( upper+lower + (upper-lower)%2 )/2;
                       }
                }
            }

        return interval;
       }

int XINTERVAL_MULTIPLESETS::XInterval(int i,double X) const
      { if(X<XMin(i) || X>XMax(i)){ throw OUT_OF_RANGE<double>(X,XMin(i),XMax(i),string("XInterval"));}

        if(interval[i-1]<1 || interval[i-1]>static_cast<int>(NX(i))-1){ interval[i-1]=(static_cast<int>(NX(i))-1)/2;}

        if(interval[i-1]!=static_cast<int>(NX(i))-1 && X>=_pX(i,interval[i-1]) && _pX(i,interval[i-1]+1)>X){ return interval[i-1];}
        if(interval[i-1]==static_cast<int>(NX(i))-1 && X>=_pX(i,interval[i-1]) && _pX(i,interval[i-1]+1)>=X){ return interval[i-1];}

        if(interval[i-1]>=2){ if(X>=_pX(i,interval[i-1]-1) && _pX(i,interval[i-1])>X){ --interval[i-1]; return interval[i-1];} }
        if(interval[i-1]<=static_cast<int>(NX(i))-3){ if(X>=_pX(i,interval[i-1]+1) && _pX(i,interval[i-1]+2)>X){ ++interval[i-1]; return interval[i-1];} }
        if(interval[i-1]==static_cast<int>(NX(i))-2){ if(X>=_pX(i,interval[i-1]+1) && _pX(i,interval[i-1]+2)>=X){ ++interval[i-1]; return interval[i-1];} }
  
        if(Equality(X,_pX(i,1))==true){ return interval[i-1]=1;} if(Equality(X,_pX(i,NX(i)))==true){ return interval[i-1]=static_cast<int>(NX(i))-1;}

        int lower=1, upper=static_cast<int>(NX(i))-1;
        interval[i-1]=( upper+lower + (upper-lower)%2 )/2;

        while(lower!=upper && (X<_pX(i,interval[i-1]) || X>=_pX(i,interval[i-1]+1)) )
             { if(X>=_pX(i,interval[i-1]+1)){ lower=interval[i-1]+1;} else{ upper=interval[i-1]-1;}
               interval[i-1]=( upper+lower + (upper-lower)%2 )/2;
              }

        return interval[i-1];
       }

// *******************************************************************

int YINTERVAL_SINGLESET::YInterval(double Y) const
    { if(Y<YMin() || Y>YMax()){ throw OUT_OF_RANGE<double>(Y,YMin(),YMax(),string("YInterval"));}

      if(interval<1 || interval>static_cast<int>(NY())-1){ interval=(static_cast<int>(NY())-1)/2;}

      if(interval!=static_cast<int>(NY())-1 && Y>=_pY(interval) && _pY(interval+1)>Y){ return interval;}
      if(interval==static_cast<int>(NY())-1 && Y>=_pY(interval) && _pY(interval+1)>=Y){ return interval;}

      if(interval>=2){ if(Y>=_pY(interval-1) && _pY(interval)>Y){ --interval; return interval;} }
      if(interval<=static_cast<int>(NY())-3){ if(Y>=_pY(interval+1) && _pY(interval+2)>Y){ ++interval; return interval;} }
      if(interval==static_cast<int>(NY())-2){ if(Y>=_pY(interval+1) && _pY(interval+2)>=Y){ ++interval; return interval;} }

      if(Equality(Y,_pY(1))==true){ return interval=1;} if(Equality(Y,_pY(NY()))==true){ return interval=static_cast<int>(NY())-1;}

      int lower=1, upper=static_cast<int>(NY())-1;
      interval=( upper+lower + (upper-lower)%2 )/2;

      while(lower!=upper && (Y<_pY(interval) || Y>=_pY(interval+1)) )
           { if(Y>=_pY(interval+1)){ lower=interval+1;} else{ upper=interval-1;}
        interval=( upper+lower + (upper-lower)%2 )/2;
       }

      return interval; 
     }

// *******************************************************************

std::vector<int> YINTERVAL_MULTIPLESETS::YInterval(std::vector<double> Y) const
      { bool found;
        for(int i=1;i<=static_cast<int>(NYSets());i++)
           { found=false;
             if(Y[i-1]<YMin(i) || Y[i-1]>YMax(i)){ throw OUT_OF_RANGE<double>(Y[i-1],YMin(i),YMax(i),string("YInterval"));}

             if(interval[i-1]<1 || interval[i-1]>static_cast<int>(NY(i))-1){ interval[i-1]=(static_cast<int>(NY(i))-1)/2;}

             if(interval[i-1]>=2){ if(Y[i-1]>=_pY(i,interval[i-1]-1) && _pY(i,interval[i-1])>Y[i-1]){ --interval[i-1]; found=true;} }
             if(interval[i-1]<=static_cast<int>(NY(i))-3){ if(Y[i-1]>=_pY(i,interval[i-1]+1) && _pY(i,interval[i-1]+2)>Y[i-1]){ ++interval[i-1]; found=true;} }
             if(interval[i-1]==static_cast<int>(NY(i))-2){ if(Y[i-1]>=_pY(i,interval[i-1]+1) && _pY(i,interval[i-1]+2)>=Y[i-1]){ ++interval[i-1]; found=true;} }
  
             if(Equality(Y[i-1],_pY(i,1))==true){ interval[i-1]=1; found=true;} 
             if(Equality(Y[i-1],_pY(i,NY(i)))==true){ interval[i-1]=static_cast<int>(NY(i))-1; found=true;}

             if(found==false)
               { int lower=1, upper=static_cast<int>(NY(i))-1;
                 interval[i-1]=( upper+lower + (upper-lower)%2 )/2;

                 while(lower!=upper && (Y[i-1]<_pY(i,interval[i-1]) || Y[i-1]>=_pY(i,interval[i-1]+1)) )
                      { if(Y[i-1]>=_pY(i,interval[i-1]+1)){ lower=interval[i-1]+1;} else{ upper=interval[i-1]-1;}
                        interval[i-1]=( upper+lower + (upper-lower)%2 )/2;
                       }
                }
            }

        return interval;
       }

int YINTERVAL_MULTIPLESETS::YInterval(int i,double Y) const
      { if(Y<YMin(i) || Y>YMax(i)){ throw OUT_OF_RANGE<double>(Y,YMin(i),YMax(i),string("YInterval"));}

        if(interval[i-1]<1 || interval[i-1]>static_cast<int>(NY(i))-1){ interval[i-1]=(static_cast<int>(NY(i))-1)/2;}

        if(interval[i-1]!=static_cast<int>(NY(i))-1 && Y>=_pY(i,interval[i-1]) && _pY(i,interval[i-1]+1)>Y){ return interval[i-1];}
        if(interval[i-1]==static_cast<int>(NY(i))-1 && Y>=_pY(i,interval[i-1]) && _pY(i,interval[i-1]+1)>=Y){ return interval[i-1];}

        if(interval[i-1]>=2){ if(Y>=_pY(i,interval[i-1]-1) && _pY(i,interval[i-1])>Y){ --interval[i-1]; return interval[i-1];} }
        if(interval[i-1]<=static_cast<int>(NY(i))-3){ if(Y>=_pY(i,interval[i-1]+1) && _pY(i,interval[i-1]+2)>Y){ ++interval[i-1]; return interval[i-1];} }
        if(interval[i-1]==static_cast<int>(NY(i))-2){ if(Y>=_pY(i,interval[i-1]+1) && _pY(i,interval[i-1]+2)>=Y){ ++interval[i-1]; return interval[i-1];} }
  
        if(Equality(Y,_pY(i,1))==true){ return interval[i-1]=1;} if(Equality(Y,_pY(i,NY(i)))==true){ return interval[i-1]=static_cast<int>(NY(i))-1;}

        int lower=1, upper=static_cast<int>(NY(i))-1;
        interval[i-1]=( upper+lower + (upper-lower)%2 )/2;

        while(lower!=upper && (Y<_pY(i,interval[i-1]) || Y>=_pY(i,interval[i-1]+1)) )
             { if(Y>=_pY(i,interval[i-1]+1)){ lower=interval[i-1]+1;} else{ upper=interval[i-1]-1;}
               interval[i-1]=( upper+lower + (upper-lower)%2 )/2;
              }

        return interval[i-1];
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

       for(i=1;i<=static_cast<int>(N())-1;i++)
          { if(Equality(_pX(i),_pX(i+1))==false){ newx.push_back(_pX(i)); newy.push_back(_pY(i));}
            else{ cout<<"\nReplacing the point "<<i<<" with x= "<<_pX(i)<<" y= "<<_pY(i);
                  int j=1; double mean=_pY(i); bool end=false;
                  while(end==false){ cout<<"\nand the point "<<i+j<<" with x= "<<_pX(i+j)<<" y= "<<_pY(i+j);
                           mean=( mean*j + _pY(i+j) )/(1.+j);
                           if(i+j>=static_cast<int>(N())){ end=true;} else{ if(Equality(_pX(i),_pX(i+j+1))==false){ end=true;} else{ j++;} }
                          };
                  cout<<"\nwith a point at x= "<<_pX(i)<<" and y= "<<mean; cout.flush();
                  newx.push_back(_pX(i)); newy.push_back(mean);
                  i+=j;
                 }
           }
       if(i==static_cast<int>(N())){ newx.push_back(_pX(N())); newy.push_back(_pY(N()));}
     
       _pSetX(newx);
       _pSetY(newy); 
      }

// *******************************************************
// *******************************************************
// *******************************************************

OPERATOR& OPERATOR::operator+=(const double &A)
            { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)+A );} 
              return (*this);
             }

OPERATOR& OPERATOR::operator-=(const double &A)
             { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)-A );} 
               return (*this);
             }

OPERATOR& OPERATOR::operator*=(const double &A)
             { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)*A );} 
               return (*this);
              }

OPERATOR& OPERATOR::operator/=(const double &A)
       { if(Equality(A,0.)==true){ throw DIVISION_BY_ZERO("/=double");}
         for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY(i,_pY(i)/A);}
         return (*this);
        }

// *******************************************************************

OPERATOR& OPERATOR::operator=(const NULLARYFUNCTOR<double> &NF)
            { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,NF() );} 
              return (*this);
             }

OPERATOR& OPERATOR::operator+=(const NULLARYFUNCTOR<double> &NF)
            { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)+NF() );}
              return (*this);
             }

OPERATOR& OPERATOR::operator-=(const NULLARYFUNCTOR<double> &NF)
             { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)-NF() );} 
               return (*this);
             }

OPERATOR& OPERATOR::operator*=(const NULLARYFUNCTOR<double> &NF)
             { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)*NF() );} 
               return (*this);
              }

OPERATOR& OPERATOR::operator/=(const NULLARYFUNCTOR<double> &NF)
       { double D;
         for(int i=1;i<=static_cast<int>(NY());i++)
            { D=NF();
               if(Equality(D,0.)==true){ throw DIVISION_BY_ZERO("/=NULLARYFUNCTOR");}
              _pSetY(i,_pY(i)/D);
             }
         return (*this);
        }

// *******************************************************************

OPERATOR& OPERATOR::operator=(const UNARYFUNCTOR<double> &UF)
             { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,UF(_pX(i)) );} 
               return (*this);
              }

OPERATOR& OPERATOR::operator+=(const UNARYFUNCTOR<double> &UF)
             { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)+UF(_pX(i)) );} 
               return (*this);
             }

OPERATOR& OPERATOR::operator-=(const UNARYFUNCTOR<double> &UF)
             { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)-UF(_pX(i)) );} 
               return (*this);
              }

OPERATOR& OPERATOR::operator*=(const UNARYFUNCTOR<double> &UF)
            { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY( i,_pY(i)*UF(_pX(i)) );} 
              return (*this);
             }

OPERATOR& OPERATOR::operator/=(const UNARYFUNCTOR<double> &UF)
        { double D;
          for(int i=1;i<=static_cast<int>(NY());i++)
                 { D=UF(_pX(i));
                     if(Equality(D,0.)==true){ throw DIVISION_BY_ZERO("/=UNARYFUNCTOR");}
                      _pSetY(i,_pY(i)/D);
                   }
          return (*this);
        }

// *******************************************************************

OPERATOR& OPERATOR::operator+=(const BINARYFUNCTOR<double> &BF)
            { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY(i,_pY(i)+BF(_pX(i),_pY(i)));} 
              return (*this);
             }

OPERATOR& OPERATOR::operator-=(const BINARYFUNCTOR<double> &BF)
            { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY(i,_pY(i)-BF(_pX(i),_pY(i)));} 
              return (*this);
             }

OPERATOR& OPERATOR::operator*=(const BINARYFUNCTOR<double> &BF)
            { for(int i=1;i<=static_cast<int>(NY());i++){ _pSetY(i,_pY(i)*BF(_pX(i),_pY(i)));} 
              return (*this);
             }

OPERATOR& OPERATOR::operator/=(const BINARYFUNCTOR<double> &BF)
       { double D;
         for(int i=1;i<=static_cast<int>(NY());i++)
                 { D=BF(_pX(i),_pY(i));
                    if(Equality(D,0.)==true){ throw DIVISION_BY_ZERO("/=BINARYFUNCTOR");}
                    _pSetY(i,_pY(i)/D);
                   }
         return (*this);
        }

} //end of namespace
