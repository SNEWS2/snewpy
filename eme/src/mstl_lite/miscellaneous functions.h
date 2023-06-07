
#include <ctime>
#include <limits>
#include <cmath>
#include <limits>
#include <vector>

#include "mstl.h"

#if !defined(_EXTRA)
#define _EXTRA

// *****************************************************************************************
// *****************************************************************************************
// *****************************************************************************************

void hold(int s);

template <typename Type> void DoNothing(Type);
template <typename Type1,typename Type2> void DoNothing(Type1,Type2);

bool Odd(int);
bool Even(int);

template <typename Type> bool Signbit(Type);
template <typename Type> Type Sign(Type);

template <typename Type> Type Plus(Type);
template <typename Type> Type Minus(Type);

template <typename Type> Type Equal(Type);
template <typename Type> Type Squared(Type);
template <typename Type> Type Cubed(Type);

template <typename Type> Type EqualDerivative(Type);
template <typename Type> Type SquaredDerivative(Type);
template <typename Type> Type CubedDerivative(Type);

template <typename Type> Type& Assign(Type&,Type);
template <typename Type> Type Add(Type,Type);
template <typename Type> Type Subtract(Type,Type);
template <typename Type> Type Multiply(Type,Type);
template <typename Type> Type Divide(Type,Type);

template <typename Type> Type Zero(Type);
template <typename Type> Type Zero(void);

template <typename Type> Type One(Type);
template <typename Type> Type One(void);

template <typename Type> Type Two(Type);
template <typename Type> Type Two(void);

#endif
