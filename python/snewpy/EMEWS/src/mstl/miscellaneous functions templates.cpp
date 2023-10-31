#if !defined(_EXTRA_FUNCTIONS)
#define _EXTRA_FUNCTIONS

// By far the stupidest function I've ever written. Still, it has its uses
template <typename Type> void DoNothing(Type T){;}
template <typename Type1,typename Type2> void DoNothing(Type1,Type2){;}

template <typename Type> bool Signbit(Type T)
{ return *(reinterpret_cast<char*>(&T)+sizeof(T)-1) & (char)128;}

template <typename Type> Type Sign(Type T){ if(Signbit(T)==true){ return -One<Type>();} return One<Type>();}

template <typename Type> Type Plus(Type T);
template <typename Type> Type Minus(Type T);

// This has to be the second stupidest function I've ever written.
template <typename Type> Type Equal(Type T){ return T;}
template <typename Type> Type Squared(Type T){ return T*T;}
template <typename Type> Type Cubed(Type T){ return T*T*T;}

template <typename Type> Type EqualDerivative(Type T){ return One(T);}
template <typename Type> Type SquaredDerivative(Type T){ return T+T;}
template <typename Type> Type CubedDerivative(Type T){ return T*T +T*T +T*T;}

template <typename Type> Type& Assign(Type &T1,Type T2){ return T1=T2;}
template <typename Type> Type Add(Type T1,Type T2){ return T1+T2;}
template <typename Type> Type Subtract(Type T1,Type T2){ return T1-T2;}
template <typename Type> Type Multiply(Type T1,Type T2){ return T1*T2;}
template <typename Type> Type Divide(Type T1,Type T2){ return T1/T2;}

template <typename Type> Type Zero(Type){ return static_cast<Type>(0);}
template <typename Type> Type Zero(void){ return static_cast<Type>(0);}

template <typename Type> std::vector<Type> Zero(std::vector<Type> &v){ return std::vector<Type>(v.size(),Zero<Type>());}

template <typename Type> std::vector<std::vector<Type> > Zero(std::vector<std::vector<Type> > &v)
         { std::vector<std::vector<Type> > w(v.size());
           for(int i=0;i<=(int)v.size()-1;i++){ w[i]=std::vector<Type>(v[i].size(),Zero<Type>());}
           return w;
          }

template <typename Type> std::vector<std::vector<std::vector<Type> > > Zero(std::vector<std::vector<std::vector<Type> > > &v)
         { std::vector<std::vector<std::vector<Type> > > w(v.size());
           for(int i=0;i<=(int)v.size()-1;i++)
              { w[i]=std::vector<std::vector<Type> >(v[i].size());
                for(int j=0;j<=(int)v[i].size()-1;j++){ w[i][j]=std::vector<Type>(v[i][j].size(),Zero<Type>());} 
               }
           return w;
          }

template <typename Type> Type One(Type){ return static_cast<Type>(1);}
template <typename Type> Type One(void){ return static_cast<Type>(1);}

template <typename Type> Type Two(Type){ return static_cast<Type>(2);}
template <typename Type> Type Two(void){ return static_cast<Type>(2);}

#endif
