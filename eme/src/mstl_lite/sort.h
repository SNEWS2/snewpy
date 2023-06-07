
#include<functional>
#include<vector>

#include "mstl.h"

#if !defined(_SORT)
#define _SORT

enum ordering { ascending, descending };

// test to see if the list is already ordered
template <template<typename,typename> class Container,typename Type,typename TypeAllocator,typename Compare> 
bool Ordered(const Container<Type,TypeAllocator>&,const Compare &C=std::less_equal<Type>());

template <template<typename,typename> class Container,typename Type,typename TypeAllocator=std::allocator<Type> > 
bool Ordered(const Container<Type,TypeAllocator>&,ordering o=ascending);

// test to see if the list is already ordered
template <typename Container,typename Type,typename Compare> 
bool Ordered(const Container&,Type,const Compare &C=std::less_equal<Type>());

template <typename Container,typename Type> 
bool Ordered(const Container&,Type,ordering o=ascending);

// *************

// sort single lists
template <template<typename,typename> class Container,typename Type,typename TypeAllocator=std::allocator<Type> > 
void StraightInsertionSort(const Container<Type,TypeAllocator> &TL,ordering o=ascending);

template <typename Container,typename Type> 
void StraightInsertionSort(const Container &TL,ordering o=ascending);

// *************

template <template<typename,typename> class Container,typename Type,typename TypeAllocator=std::allocator<Type> > 
void ShellSort(const Container<Type,TypeAllocator> &TL,ordering o=ascending);

template <typename Container,typename Type> 
void ShellSort(const Container &TL,ordering o=ascending);

// *************

// return a list in decreasing order of the magnitude of the elements, I being the order of the points in the original list
template <template<typename,typename> class TContainer,class IContainer,typename Type,typename TypeAllocator> 
TContainer<Type,TypeAllocator> LejaSort(TContainer<Type,TypeAllocator> L1,IContainer &I);

template <typename TContainer,typename IContainer,typename Type> 
TContainer LejaSort(TContainer L1,IContainer &I,Type);

// *********************************************************************************************************
// *********************************************************************************************************
// *********************************************************************************************************

// sort two lists, the first is the key
template <template<typename,typename> class Container1,template<typename,typename> class Container2,typename Type1,typename Type2,typename Type1Allocator,typename Type2Allocator,typename Compare> 
void Sort(Container1<Type1,Type1Allocator> &c1,Container2<Type2,Type2Allocator> &c2,const Compare &C=std::less_equal<Type1>());

template <template<typename,typename> class Container1,template<typename,typename> class Container2,typename Type1,typename Type2=Type1,typename Type1Allocator=std::allocator<Type1>,typename Type2Allocator=std::allocator<Type2> > 
void Sort(Container1<Type1,Type1Allocator> &c1,Container2<Type2,Type2Allocator> &c2,ordering o);

// *************

template <typename Container1,typename Container2,typename Type1,typename Type2,typename Compare> 
void Sort(Container1 &c1,Container2 &c2,Type1 t1,Type2 t2,const Compare &C=std::less_equal<Type1>());

template <typename Container1,typename Container2,typename Type1,typename Type2> 
void Sort(Container1 &c1,Container2 &c2,Type1 t1,Type2 t2,ordering o);

// *********************************************************************************************************
// *********************************************************************************************************
// *********************************************************************************************************

template <template<typename,typename> class Container1,template<typename,typename> class Container2,typename Type1,typename Type2,typename Type1Allocator,typename Type2Allocator,typename Compare>
void StraightInsertionSort(Container1<Type1,Type1Allocator> &c1,Container2<Type2,Type2Allocator> &c2,const Compare &C=std::less_equal<Type1>());

template <template<typename,typename> class Container1,template<typename,typename> class Container2,typename Type1,typename Type2,typename Type1Allocator,typename Type2Allocator,typename Compare>
void StraightInsertionSort(Container1<Type1,Type1Allocator> &c1,Container2<Type2,Type2Allocator> &c2,ordering o);

// *************

template <typename Container1,typename Container2,typename Type1,typename Type2,typename Compare>
void StraightInsertionSort(Container1 &c1,Container2 &c2,Type1,Type2,const Compare &C);

template <typename Container1,typename Container2,typename Type1,typename Type2>
void StraightInsertionSort(Container1 &c1,Container2 &c2,Type1 t1,Type2 t2,ordering o);

#endif
