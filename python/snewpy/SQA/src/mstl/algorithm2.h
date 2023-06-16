
#include <algorithm>
#include <iterator>

#include "mstl.h"

#if !defined(_ALGORITHM2)
#define _ALGORITHM2

template <class InputIterator,class OutputIterator>
OutputIterator catenate(InputIterator first1,InputIterator last1,InputIterator first2,InputIterator last2,OutputIterator result);

#endif
