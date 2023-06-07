
template <class InputIterator,class OutputIterator>
OutputIterator catenate(InputIterator first1,InputIterator last1,InputIterator first2,InputIterator last2,OutputIterator result)
         { return copy(first2,last2,copy(first1,last1,result));}

