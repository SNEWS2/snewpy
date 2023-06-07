#if !defined(_SORT_FUNCTIONS)
#define _SORT_FUNCTIONS

template <template<typename,typename> class Container,typename Type,typename TypeAllocator,typename Compare>
bool Ordered(const Container<Type,TypeAllocator> &c,const Compare &C)
     { bool ordered=true;
       int i=0;
       while(ordered==true && i<=(int)c.size()-2)
            { if( C(c[i],c[i+1])==false ){ ordered=false;} i++;}
       return ordered;
      }

template <template<typename,typename> class Container,typename Type,typename TypeAllocator>
bool Ordered(const Container<Type,TypeAllocator> &c,ordering o)
     { if(o==ascending){ return Ordered(c,std::less_equal<Type>());}
       if(o==descending){ return Ordered(c,std::greater_equal<Type>());}
       throw UNKNOWN_ERROR("Ordered(Container,ordering)");
      }

template <typename Container,typename Type,typename Compare>
bool Ordered(const Container &c,Type t,const Compare &C)
     { bool ordered=true;
       int i=0;
       while(ordered==true && i<=(int)c.size()-2)
            { if( C(c[i],c[i+1])==false ){ ordered=false;} i++;}
       return ordered;
      }

template <typename Container,typename Type>
bool Ordered(const Container &c,Type t,ordering o)
     { if(o==ascending){ return Ordered(c,t,std::less_equal<Type>());}
       if(o==descending){ return Ordered(c,t,std::greater_equal<Type>());}
       throw UNKNOWN_ERROR("Ordered(Container,ordering)");
      }

//******************************

template <template<typename,typename> class Container,typename Type,typename TypeAllocator> 
void StraightInsertionSort(const Container<Type,TypeAllocator> &TL,ordering o)
         { if(Ordered(TL,o)==true){ return;}
           Type T;
           for(int j=1;j<=(int)TL.size()-1;j++)
              { T=TL[j];
                int i=j-1;
                if(o==ascending){ while(i>=0 && TL[i]>T){ TL[i+1]=TL[i]; i--;}; }
                if(o==descending){ while(i>=0 && TL[i]<T){ TL[i+1]=TL[i]; i--;}; }
                TL[i+1]=T;
               }
          }

//specialization for arrays of pointers
template <template<typename,typename> class Container,typename Type,typename TypeAllocator> 
void StraightInsertionSort(int N,const Container<Type,TypeAllocator> &TL,ordering o)
         { Type *T;
           for(int j=1;j<=(int)TL.size()-1;j++)
              { (*T)=(*TL[j]);
                int i=j-1;
                if(o==ascending){ while(i>=0 && (*TL[i])>(*T)){ (*TL[i+1])=(*TL[i]); i--;}; }
                if(o==descending){ while(i>=0 && (*TL[i])<(*T)){ (*TL[i+1])=(*TL[i]); i--;}; }
                *TL[i+1]=*T;
               }
          }

template <typename Container,typename Type> 
void StraightInsertionSort(const Container &TL,Type t,ordering o)
         { if(Ordered(TL,t,o)==true){ return;}
           Type T;
           for(int j=1;j<=(int)TL.size()-1;j++)
              { T=TL[j];
                int i=j-1;
                if(o==ascending){ while(i>=0 && TL[i]>T){ TL[i+1]=TL[i]; i--;}; }
                if(o==descending){ while(i>=0 && TL[i]<T){ TL[i+1]=TL[i]; i--;}; }
                TL[i+1]=T;
               }
          }

//specialization for arrays of pointers
template <typename Container,typename Type>
void StraightInsertionSort(int N,const Container &TL,Type *t,ordering o)
         { Type *T;
           for(int j=1;j<=(int)TL.size()-1;j++)
              { (*T)=(*TL[j]);
                int i=j-1;
                if(o==ascending){ while(i>=0 && (*TL[i])>(*T)){ (*TL[i+1])=(*TL[i]); i--;}; }
                if(o==descending){ while(i>=0 && (*TL[i])<(*T)){ (*TL[i+1])=(*TL[i]); i--;}; }
                *TL[i+1]=*T;
               }
          }

//******************************

template <template<typename,typename> class Container,typename Type,typename TypeAllocator> 
void ShellSort(const Container<Type,TypeAllocator> &TL,ordering o)
         { if(Ordered(TL,o)==true){ return;}
           Type T;
           int i,j,gap,n = static_cast<int>(TL.size());

           for(gap=n/2;gap>=1;gap/=2)
              { for(i=gap;i<=n-1;i++)
                   { T=TL[i];
                     if(o==ascending){ for(j=i;j>=gap && TL[j-gap]>T;j-=gap){ TL[j]=TL[j-gap];} }
                     if(o==descending){ for(j=i;j>=gap && TL[j-gap]<T;j-=gap){ TL[j]=TL[j-gap];} }
                     TL[j] = T;
                    }
                }
          }

template <typename Container,typename Type> 
void ShellSort(const Container &TL,Type t,ordering o)
         { if(Ordered(TL,t,o)==true){ return;}
           Type T;
           int i,j,gap,n = static_cast<int>(TL.size());

           for(gap=n/2;gap>=1;gap/=2)
              { for(i=gap;i<=n-1;i++)
                   { T=TL[i];
                     if(o==ascending){ for(j=i;j>=gap && TL[j-gap]>T;j-=gap){ TL[j]=TL[j-gap];} }
                     if(o==descending){ for(j=i;j>=gap && TL[j-gap]<T;j-=gap){ TL[j]=TL[j-gap];} }
                     TL[j] = T;
                    }
                }
          }

//************************

template <template<typename,typename> class TContainer,class IContainer,typename Type,typename TypeAllocator> 
TContainer<Type,TypeAllocator> LejaSort(TContainer<Type,TypeAllocator> L1,IContainer &I)
         { TContainer<Type,TypeAllocator> L2(L1.size());
           std::vector<bool> used(L1.size());

           L2[0]=L1[0]; I[0]=0; used[0]=false;
           for(int i=1;i<=(int)L1.size()-1;i++){ used[i]=false; if(fabs(L1[i])>fabs(L2[0])){ L2[0]=L1[i]; I[0]=i;} }
           used[I[0]]=true;

           double magj, magk;
           for(int i=1;i<=(int)L1.size()-1;i++)
              { L2[i]=L1[i]; I[i]=i;
                magk=1.; for(int k=0;k<=i-1;k++){ magk*=fabs(L2[i]-L2[k]);}

                for(int j=0;j<=(int)L1.size()-1;j++)
                   { if(used[j]==false){ magj=1.; for(int k=0;k<=i-1;k++){ magj*=fabs(L1[j]-L2[k]);}
                                         if(magj>magk){ L2[i]=L1[j]; I[i]=j; magk=magj;}
                                        }
                    }
                used[I[i]]=true;
               }

           return L2;
          }

template <typename TContainer,typename IContainer,typename Type> 
TContainer LejaSort(TContainer L1,IContainer &I,Type T)
         { TContainer L2(L1.size());
           std::vector<bool> used(L1.size());

           L2[0]=L1[0]; I[0]=0; used[0]=false;
           for(int i=1;i<=(int)L1.size()-1;i++){ used[i]=false; if(fabs(L1[i])>fabs(L2[0])){ L2[0]=L1[i]; I[0]=i;} }
           used[I[0]]=true;

           double magj, magk;
           for(int i=1;i<=(int)L1.size()-1;i++)
              { L2[i]=L1[i]; I[i]=i;
                magk=1.; for(int k=0;k<=i-1;k++){ magk*=fabs(L2[i]-L2[k]);}

                for(int j=0;j<=(int)L1.size()-1;j++)
                   { if(used[j]==false){ magj=1.; for(int k=0;k<=i-1;k++){ magj*=fabs(L1[j]-L2[k]);}
                                         if(magj>magk){ L2[i]=L1[j]; I[i]=j; magk=magj;}
                                        }
                    }
                used[I[i]]=true;
               }

           return L2;
          }

//*****************************************************************************************
//*****************************************************************************************
//*****************************************************************************************

// Two lists
template <template<typename,typename> class Container1,template<typename,typename> class Container2,typename Type1,typename Type2,typename Type1Allocator,typename Type2Allocator,typename Compare>
void Sort(Container1<Type1,Type1Allocator> &c1,Container2<Type2,Type2Allocator> &c2,const Compare &C)
         { try{ StraightInsertionSort(c1,c2,C);}
           catch(DIFFERENT_LENGTHS &DL){ DL.ChangeFunction("Sort(Container,Container,Compare)"); throw DL;}
          }

template <template<typename,typename> class Container1,template<typename,typename> class Container2,typename Type1,typename Type2,typename Type1Allocator,typename Type2Allocator>
void Sort(Container1<Type1,Type1Allocator> &c1,Container2<Type2,Type2Allocator> &c2,ordering o)
         { try{ if(o==ascending){ Sort(c1,c2,std::less_equal<Type1>()); return; }
                if(o==descending){ Sort(c1,c2,std::greater_equal<Type1>()); return;}
                throw UNKNOWN_ERROR("Sort(Container,Container,ordering)");
               }
           catch(DIFFERENT_LENGTHS &DL){ DL.ChangeFunction("Sort(Container,Container,ordering)"); throw DL;}
          }

//********************************

template <typename Container1,typename Container2,typename Type1,typename Type2,typename Compare>
void Sort(Container1 &c1,Container2 &c2,Type1 t1,Type2 t2,const Compare &C)
         { try{ StraightInsertionSort(c1,c2,t1,t2,C);}
           catch(DIFFERENT_LENGTHS &DL){ DL.ChangeFunction("Sort(Container,Container,Compare)"); throw DL;}
          }

template <typename Container1,typename Container2,typename Type1,typename Type2>
void Sort(Container1 &c1,Container2 &c2,Type1 t1,Type2 t2,ordering o)
         { try{ if(o==ascending){ Sort(c1,c2,t1,t2,std::less_equal<Type1>()); return; }
                if(o==descending){ Sort(c1,c2,t1,t2,std::greater_equal<Type1>()); return;}
                throw UNKNOWN_ERROR("Sort(Container,Container,ordering)");
               }
           catch(DIFFERENT_LENGTHS &DL){ DL.ChangeFunction("Sort(Container,Container,ordering)"); throw DL;}
          }

//*****************************************************************************************
//*****************************************************************************************
//*****************************************************************************************

template <template<typename,typename> class Container1,template<typename,typename> class Container2,typename Type1,typename Type2,typename Type1Allocator,typename Type2Allocator,typename Compare>
void StraightInsertionSort(Container1<Type1,Type1Allocator> &c1,Container2<Type2,Type2Allocator> &c2,const Compare &C)
         { if(c1.size()!=c2.size()){ throw DIFFERENT_LENGTHS("StraightInsertionSort(Container,Container,Compare)");}
           if(Ordered(c1,C)==true){ return;}
           Type1 T1; Type2 T2;
           for(int j=1;j<=(int)c1.size()-1;j++)
              { T1=c1[j]; T2=c2[j];
                int i=j-1;
                while (i>=0 && C(T1,c1[i])==true){ c1[i+1]=c1[i]; c2[i+1]=c2[i]; i--;};
                c1[i+1]=T1; c2[i+1]=T2;
               }
          }

template <template<typename,typename> class Container1,template<typename,typename> class Container2,typename Type1,typename Type2,typename Type1Allocator,typename Type2Allocator>
void StraightInsertionSort(Container1<Type1,Type1Allocator> &c1,Container2<Type2,Type2Allocator> &c2,ordering o)
         { try{ if(o==ascending){ StraightInsertionSort(c1,c2,std::less_equal<Type1>()); return; }
                if(o==descending){ StraightInsertionSort(c1,c2,std::greater_equal<Type1>()); return;}
                throw UNKNOWN_ERROR("StraightInsertionSort(Container,Container,ordering)");
               }
           catch(DIFFERENT_LENGTHS &DL){ DL.ChangeFunction("StraightInsertionSort(Container,Container,ordering)"); throw DL;}
          }

//*****************

template <typename Container1,typename Container2,typename Type1,typename Type2,typename Compare>
void StraightInsertionSort(Container1 &c1,Container2 &c2,Type1 t1,Type2 t2,const Compare &C)
         { if(c1.size()!=c2.size()){ throw DIFFERENT_LENGTHS("StraightInsertionSort(Container,Container,Compare)");}
           if(Ordered(c1,t1,C)==true){ return;}
           Type1 T1; Type2 T2;
           for(int j=1;j<=(int)c1.size()-1;j++)
              { T1=c1[j]; T2=c2[j];
                int i=j-1;
                while (i>=0 && C(T1,c1[i])==true){ c1[i+1]=c1[i]; c2[i+1]=c2[i]; i--;};
                c1[i+1]=T1; c2[i+1]=T2;
               }
          }

template <typename Container1,typename Container2,typename Type1,typename Type2>
void StraightInsertionSort(Container1 &c1,Container2 &c2,Type1 t1,Type2 t2,ordering o)
         { try{ if(o==ascending){ StraightInsertionSort(c1,c2,t1,t2,std::less_equal<Type1>()); return; }
                if(o==descending){ StraightInsertionSort(c1,c2,t1,t2,std::greater_equal<Type1>()); return;}
                throw UNKNOWN_ERROR("StraightInsertionSort(Container,Container,ordering)");
               }
           catch(DIFFERENT_LENGTHS &DL){ DL.ChangeFunction("StraightInsertionSort(Container,Container,ordering)"); throw DL;}
          }

#endif
