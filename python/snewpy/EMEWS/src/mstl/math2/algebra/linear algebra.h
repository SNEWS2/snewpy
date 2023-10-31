
#include<limits>
#include<cmath>

#include "mstl.h"

// *******************************************************************
// *******************************************************************
// *******************************************************************

#if !defined(_LINALG)
#define _LINALG

template <typename Type,typename VType> struct EIGENPAIR;

// *******************************************************************
// forward declarations

template <typename Type,typename MType> class MATRIXEXPRESSION;
template <typename Type,std::size_t N1,std::size_t N2> class MATRIX;

template <typename Type,typename CVType> class CVECTOREXPRESSION;
template <typename Type,std::size_t N> class CVECTOR;

template <typename Type,typename RVType> class RVECTOREXPRESSION;

template <typename Type,std::size_t N> class RVECTOR;

template <typename Type,typename MVType> class MVECTOREXPRESSION;
template <typename Type,std::size_t N> class MVECTOR;

// *******************************************************************
// *******************************************************************
// ************************* FUNCTION PROTYPES ***********************
// *******************************************************************
// *******************************************************************

template <typename Type,typename MType> MATRIXEXPRESSION_NEGATE<Type,MATRIXEXPRESSION<Type,MType> > operator-(MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIXEXPRESSION_MATRIXTRANSPOSE<Type,MATRIXEXPRESSION<Type,MType> > Transpose(MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIXEXPRESSION_MATRIXANTITRANSPOSE<Type,MATRIXEXPRESSION<Type,MType> > AntiTranspose(MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIXEXPRESSION_MATRIXADJOINT<Type,MATRIXEXPRESSION<Type,MType> > Adjoint(MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIXEXPRESSION_CONJUGATE<Type,MATRIXEXPRESSION<Type,MType> > Conjugate(MATRIXEXPRESSION<Type,MType> const&);

// *******************************************************************

template <typename Type,typename MType1,typename MType2> MATRIXEXPRESSION_MATRIXADDITION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> > operator+(MATRIXEXPRESSION<Type,MType1> const&,MATRIXEXPRESSION<Type,MType2> const&);

template <typename Type,typename MType1,typename MType2> MATRIXEXPRESSION_MATRIXSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> > operator-(MATRIXEXPRESSION<Type,MType1> const&,MATRIXEXPRESSION<Type,MType2> const&);

template <typename Type,typename MType1,typename MType2> MATRIXEXPRESSION_MATRIXMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> > operator*(MATRIXEXPRESSION<Type,MType1> const&,MATRIXEXPRESSION<Type,MType2> const&);

template <typename Type,typename MType1,typename MType2> MATRIXEXPRESSION_MATRIXMULTIPLICATION<std::complex<Type>,MATRIXEXPRESSION<std::complex<Type>,MType1>,MATRIXEXPRESSION<Type,MType2> > operator*(MATRIXEXPRESSION<std::complex<Type>,MType1> const&,MATRIXEXPRESSION<Type,MType2> const&);

template <typename Type,typename MType1,typename MType2> MATRIXEXPRESSION_MATRIXMULTIPLICATION<std::complex<Type>,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<std::complex<Type>,MType2> > operator*(MATRIXEXPRESSION<Type,MType1> const&,MATRIXEXPRESSION<std::complex<Type>,MType2> const&);

//template <typename Type,typename MType1,typename MType2> MATRIXEXPRESSION_MATRIXDIVISION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> > operator/(MATRIXEXPRESSION<Type,MType1> const&,MATRIXEXPRESSION<Type,MType2> const&);

// *******************************************************************

template <typename Type,typename MType> MATRIXEXPRESSION_SCALARMATRIXADDITION<Type,MATRIXEXPRESSION<Type,MType> > operator+(Type const&,MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIXEXPRESSION_MATRIXSCALARADDITION<Type,MATRIXEXPRESSION<Type,MType> > operator+(MATRIXEXPRESSION<Type,MType> const&,Type const&);

template <typename Type,typename MType> MATRIXEXPRESSION_SCALARMATRIXSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType> > operator-(Type const&,MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIXEXPRESSION_MATRIXSCALARSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType> > operator-(MATRIXEXPRESSION<Type,MType> const&,Type const&);

template <typename Type,typename MType> MATRIXEXPRESSION_SCALARMATRIXMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType> > operator*(Type const&,MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIXEXPRESSION_MATRIXSCALARMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType> > operator*(MATRIXEXPRESSION<Type,MType> const&,Type const&);

//template <typename Type,typename MType> MATRIXEXPRESSION_SCALARMATRIXDIVISION<Type,MATRIXEXPRESSION<Type,MType> > operator/(Type const&,MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIXEXPRESSION_MATRIXSCALARDIVISION<Type,MATRIXEXPRESSION<Type,MType> > operator/(MATRIXEXPRESSION<Type,MType> const&,Type const&);

// *******************************************************************

template <typename Type,typename CVType,typename RVType> MATRIXEXPRESSION_VECTORPRODUCT<Type,CVECTOREXPRESSION<Type,CVType>,RVECTOREXPRESSION<Type,RVType> > operator*(CVECTOREXPRESSION<Type,CVType> const&,RVECTOREXPRESSION<Type,RVType> const&);

// *******************************************************************
// *******************************************************************
// *******************************************************************
// *******************************************************************
// *******************************************************************

template <typename Type,typename CVType> CVECTOREXPRESSION_NEGATE<Type,CVECTOREXPRESSION<Type,CVType> > operator-(CVECTOREXPRESSION<Type,CVType> const&);

template <typename Type,typename CVType> RVECTOREXPRESSION_CVECTORTRANSPOSE<Type,CVECTOREXPRESSION<Type,CVType> > Transpose(CVECTOREXPRESSION<Type,CVType> const&);

template <typename Type,typename CVType> RVECTOREXPRESSION_CVECTORADJOINT<Type,CVECTOREXPRESSION<Type,CVType> > Adjoint(CVECTOREXPRESSION<Type,CVType> const&);

template <typename Type,typename CVType> CVECTOREXPRESSION_CONJUGATE<Type,CVECTOREXPRESSION<Type,CVType> > Conjugate(CVECTOREXPRESSION<Type,CVType> const&);

// *******************************************************************

template <typename Type,typename CVType1,typename CVType2> CVECTOREXPRESSION_CVECTORADDITION<Type,CVECTOREXPRESSION<Type,CVType1>,CVECTOREXPRESSION<Type,CVType2> > operator+(CVECTOREXPRESSION<Type,CVType1> const&,CVECTOREXPRESSION<Type,CVType2> const&);

template <typename Type,typename CVType1,typename CVType2> CVECTOREXPRESSION_CVECTORSUBTRACTION<Type,CVECTOREXPRESSION<Type,CVType1>,CVECTOREXPRESSION<Type,CVType2> > operator-(CVECTOREXPRESSION<Type,CVType1> const&,CVECTOREXPRESSION<Type,CVType2> const&);

// *******************************************************************

template <typename Type,typename CVType> CVECTOREXPRESSION_CVECTORSCALARMULTIPLICATION<Type,CVECTOREXPRESSION<Type,CVType> > operator*(CVECTOREXPRESSION<Type,CVType> const&,Type const&);

template <typename Type,typename CVType> CVECTOREXPRESSION_SCALARCVECTORMULTIPLICATION<Type,CVECTOREXPRESSION<Type,CVType> > operator*(Type const&,CVECTOREXPRESSION<Type,CVType> const&);

template <typename Type,typename CVType> CVECTOREXPRESSION_CVECTORSCALARDIVISION<Type,CVECTOREXPRESSION<Type,CVType> > operator/(CVECTOREXPRESSION<Type,CVType> const&,Type const&);

// *******************************************************************

template <typename Type,typename MType,typename CVType> CVECTOREXPRESSION_MATRIXCVECTORMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType> ,CVECTOREXPRESSION<Type,CVType> > operator*(MATRIXEXPRESSION<Type,MType> const&,CVECTOREXPRESSION<Type,CVType> const&);

// *******************************************************************
// *******************************************************************
// *******************************************************************
// *******************************************************************
// *******************************************************************

template <typename Type,typename RVType> RVECTOREXPRESSION_NEGATE<Type,RVECTOREXPRESSION<Type,RVType> > operator-(RVECTOREXPRESSION<Type,RVType> const&);

template <typename Type,typename RVType> CVECTOREXPRESSION_RVECTORTRANSPOSE<Type,RVECTOREXPRESSION<Type,RVType> > Transpose(RVECTOREXPRESSION<Type,RVType> const&);

template <typename Type,typename RVType> CVECTOREXPRESSION_RVECTORADJOINT<Type,RVECTOREXPRESSION<Type,RVType> > Adjoint(RVECTOREXPRESSION<Type,RVType> const&);

template <typename Type,typename RVType> RVECTOREXPRESSION_CONJUGATE<Type,RVECTOREXPRESSION<Type,RVType> > Conjugate(RVECTOREXPRESSION<Type,RVType> const&);

// *******************************************************************

template <typename Type,typename RVType1,typename RVType2> RVECTOREXPRESSION_RVECTORADDITION<Type,RVECTOREXPRESSION<Type,RVType1>,RVECTOREXPRESSION<Type,RVType2> > operator+(RVECTOREXPRESSION<Type,RVType1> const&,RVECTOREXPRESSION<Type,RVType2> const&);

template <typename Type,typename RVType1,typename RVType2> RVECTOREXPRESSION_RVECTORSUBTRACTION<Type,RVECTOREXPRESSION<Type,RVType1>,RVECTOREXPRESSION<Type,RVType2> > operator-(RVECTOREXPRESSION<Type,RVType1> const&,RVECTOREXPRESSION<Type,RVType2> const&);

// *******************************************************************

template <typename Type,typename RVType> RVECTOREXPRESSION_RVECTORSCALARMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType> > operator*(RVECTOREXPRESSION<Type,RVType> const&,Type const&);

template <typename Type,typename RVType> RVECTOREXPRESSION_SCALARRVECTORMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType> > operator*(Type const&,RVECTOREXPRESSION<Type,RVType> const&);

template <typename Type,typename RVType> RVECTOREXPRESSION_RVECTORSCALARDIVISION<Type,RVECTOREXPRESSION<Type,RVType> > operator/(RVECTOREXPRESSION<Type,RVType> const&,Type const&);

// *******************************************************************

template <typename Type,typename RVType,typename MType> RVECTOREXPRESSION_RVECTORMATRIXMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType> ,MATRIXEXPRESSION<Type,MType> > operator*(RVECTOREXPRESSION<Type,RVType> const&,MATRIXEXPRESSION<Type,MType> const&);

// *******************************************************************

//dot product
template <typename Type,typename RVType,typename CVType> 
Type operator*(RVECTOREXPRESSION<Type,RVType> const&,CVECTOREXPRESSION<Type,CVType> const&);

// *******************************************************************
// ********************* << operator *********************************
// *******************************************************************

template <typename Type,typename MType> std::ostream& operator<<(std::ostream&,MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename CVType> std::ostream& operator<<(std::ostream&,CVECTOREXPRESSION<Type,CVType> const&);
template <typename Type,typename RVType> std::ostream& operator<<(std::ostream&,RVECTOREXPRESSION<Type,RVType> const&);

// *******************************************************************
// *******************************************************************
// ****************** Submatrices ************************************
// *******************************************************************
// *******************************************************************

template <typename Type,typename MType> MATRIX<Type,0,0> SubMatrix(MATRIXEXPRESSION<Type,MType> const&,int,int,int,int);

template <typename Type,typename MType> CVECTOR<Type,0> SubColumnVector(MATRIXEXPRESSION<Type,MType> const&,int,int,int);
template <typename Type,typename CVType> CVECTOR<Type,0> SubColumnVector(CVECTOREXPRESSION<Type,CVType> const&,int,int);

template <typename Type,typename MType> RVECTOR<Type,0> SubRowVector(MATRIXEXPRESSION<Type,MType> const&,int,int,int);
template <typename Type,typename RVType> RVECTOR<Type,0> SubRowVector(RVECTOREXPRESSION<Type,RVType> const&,int,int);

template <typename Type,typename MType1,typename MType2> 
MATRIXEXPRESSION_SUBMATRIXADDITION<Type,MType1,MType2> SubMatrixAddition(MATRIXEXPRESSION<Type,MType1> const&,MATRIXEXPRESSION<Type,MType2> const&,int,int);

template <typename Type,typename MType1,typename MType2> 
MATRIXEXPRESSION_SUBMATRIXSUBTRACTION<Type,MType1,MType2> SubMatrixSubtraction(MATRIXEXPRESSION<Type,MType1> const&,MATRIXEXPRESSION<Type,MType2> const&,int,int);

template <typename Type,typename MType,typename CVType> MATRIX<Type,0,0> Deflate(MATRIXEXPRESSION<Type,MType> const&,CVECTOREXPRESSION<Type,CVType> const&);
template <typename Type,typename MType,typename RVType> MATRIX<Type,0,0> Deflate(MATRIXEXPRESSION<Type,MType> const&,RVECTOREXPRESSION<Type,RVType> const&);

// *******************************************************************
// *******************************************************************
// ********************* SPECIAL MATRICES AND VECTORS ****************
// *******************************************************************
// *******************************************************************

template <typename Type,std::size_t N1,std::size_t N2> MATRIX<Type,N1,N2> ZeroMatrix(void);
template <typename Type> MATRIX<Type,0,0> ZeroMatrix(std::size_t N1,std::size_t N2);
template <typename Type> MATRIX<Type,0,0> ZeroMatrix(std::size_t N);

template <typename Type,std::size_t N> MATRIX<Type,N,N> UnitMatrix(void);
template <typename Type> MATRIX<Type,0,0> UnitMatrix(std::size_t N);

template <typename Type,std::size_t N> MATRIX<Type,N,N> PermutationMatrix(int i,int j);
template <typename Type> MATRIX<Type,0,0> PermutationMatrix(std::size_t N,int i,int j);

template <typename Type,std::size_t N> MATRIX<Type,N,N> ProjectionMatrix(int i);
template <typename Type> MATRIX<Type,0,0> ProjectionMatrix(std::size_t N,int i);

template <typename Type,std::size_t N> MATRIX<Type,N,N> RightShiftMatrix(void);
template <typename Type> MATRIX<Type,0,0> RightShiftMatrix(std::size_t N);

template <typename Type,typename CVType> MATRIX<Type,0,0> HouseholderMatrix(CVECTOREXPRESSION<Type,CVType> &CV);
template <typename Type,typename RVType> MATRIX<Type,0,0> HouseholderMatrix(RVECTOREXPRESSION<Type,RVType> &RV);

// ***************

template <std::size_t N> MATRIX<double,N,N> Givens(int N1,int N2,double ANGLE);
MATRIX<double,0,0> Givens(std::size_t N,int N1,int N2,double ANGLE);

template <std::size_t N> MATRIX<double,N,N> Givens(int N1,int N2,double COS,double SIN);
MATRIX<double,0,0> Givens(std::size_t N,int N1,int N2,double COS,double SIN);

template <std::size_t N> MATRIX<std::complex<double>,N,N> ComplexGivens(int N1,int N2,double ANGLE,double PHASE);
MATRIX<std::complex<double>,0,0> ComplexGivens(std::size_t N,int N1,int N2,double ANGLE,double PHASE);

template <std::size_t N> MATRIX<std::complex<double>,N,N> ComplexGivens(int N1,int N2,double COS,double SIN,double PHASE);
MATRIX<std::complex<double>,0,0> ComplexGivens(std::size_t N,int N1,int N2,double COS,double SIN,double PHASE);

template <std::size_t N> MATRIX<std::complex<double>,N,N> ComplexGivens(int N1,int N2,std::complex<double> ALPHA,std::complex<double> BETA);
MATRIX<std::complex<double>,0,0> ComplexGivens(std::size_t N,int N1,int N2,std::complex<double> ALPHA,std::complex<double> BETA);

// ***************

template <std::size_t N> MATRIX<double,N,N> IntegerMatrix(void);   // a diagonal matrix with entries equal to the integers
MATRIX<double,0,0> IntegerMatrix(std::size_t N);   // a diagonal matrix with entries equal to the integers

template <std::size_t N> MATRIX<double,N,N> FactorialMatrix(void); // a diagonal matrix with entries equal to the factorial numbers
MATRIX<double,0,0> FactorialMatrix(std::size_t N); // a diagonal matrix with entries equal to the factorial numbers

template <typename Type> MATRIX<Type,0,0> DiagonalMatrix(std::vector<Type> const &X);

template <typename Type> MATRIX<Type,0,0> VandermondeMatrix(std::vector<Type> const &X);

template <typename Type,std::size_t N> MATRIX<Type,N,N> TranslationMatrix(Type const &X); // an upper triangular matrix that shifts Vandermonde matrices
template <typename Type> MATRIX<Type,0,0> TranslationMatrix(std::size_t N,Type const &X); // an upper triangular matrix that shifts Vandermonde matrices

template <typename Type,std::size_t M> CVECTOR<Type,M> BasisColumnVector(int N);
template <typename Type> CVECTOR<Type,0> BasisColumnVector(std::size_t M,int N);

template <typename Type,std::size_t M> RVECTOR<Type,M> BasisRowVector(int N);
template <typename Type> RVECTOR<Type,0> BasisRowVector(std::size_t M,int N);

// *******************************************************************
// ***************** MATRIX PROPERTY TESTS ***************************
// *******************************************************************

template <typename Type,typename MType> bool ZeroTest(MATRIXEXPRESSION<Type,MType> const&);                // returns true if every element is zero

template <typename Type,typename MType> bool DiagonalTest(MATRIXEXPRESSION<Type,MType> const&);                               // true if only diagonal elements of a sqaure matrix are non zero
template <typename Type,typename MType> bool TridiagonalTest(MATRIXEXPRESSION<Type,MType> const&);                            // return true if matrix is tridiagonal
template <typename Type,typename MType> bool BanddiagonalTest(MATRIXEXPRESSION<Type,MType> const&,int,int); // return true if matrix is tridiagonal

template <typename Type,typename MType> bool SpurTest(MATRIXEXPRESSION<Type,MType> const&);                // return true if no spur element is zero
template <typename Type,typename MType> bool SquareTest(MATRIXEXPRESSION<Type,MType> const&);              // return true if matrix is square
template <typename Type,typename MType> bool LowerTriangleTest(MATRIXEXPRESSION<Type,MType> const&);       // returns true if a square matrix has nozero elements in the lower triangle
template <typename Type,typename MType> bool UpperTriangleTest(MATRIXEXPRESSION<Type,MType> const&);       // returns true if a square matrix has nozero elements in the upper triangle
template <typename Type,typename MType> bool SymmetricTest(MATRIXEXPRESSION<Type,MType> const&);      // returns true if a square matrix is symmetric
template <typename Type,typename MType> bool AntiSymmetricTest(MATRIXEXPRESSION<Type,MType> const&);  // returns true if a square matrix is antisymmetric

template <typename Type,typename MType> bool LowerHessenbergTest(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> bool LowerHessenbergTest(MATRIXEXPRESSION<std::complex<Type>,MType> const&);
template <typename Type,typename MType> bool UpperHessenbergTest(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> bool UpperHessenbergTest(MATRIXEXPRESSION<std::complex<Type>,MType> const&);

template <typename Type,typename MType> bool RowDiagonalDominanceTest(MATRIXEXPRESSION<Type,MType> const&);  // returns true if the square of the diagonal element of every row is greater than the sum of the squares of the rest
template <typename Type,typename MType> bool RowDiagonalDominanceTest(MATRIXEXPRESSION<std::complex<Type>,MType> const&);  
template <typename Type,typename MType> bool ColumnDiagonalDominanceTest(MATRIXEXPRESSION<Type,MType> const&);  // returns true if the square of the diagonal element of every column is greater than the sum of the squares of the rest
template <typename Type,typename MType> bool ColumnDiagonalDominanceTest(MATRIXEXPRESSION<std::complex<Type>,MType> const&); 

// *******************************************************************
// *******************************************************************
// ***************** MATRIX PROPERTIES *******************************
// *******************************************************************
// *******************************************************************

template <typename Type,typename MType> Type FrobeniusNorm(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> Type FrobeniusNorm(MATRIXEXPRESSION<std::complex<Type>,MType> const&);

template <typename Type,typename MType> MATRIX<Type,0,0> Diagonal(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> LowerTriangle(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> UpperTriangle(MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> CVECTOR<Type,0> Column(MATRIXEXPRESSION<Type,MType> const&,int);
template <typename Type,typename MType> RVECTOR<Type,0> Row(MATRIXEXPRESSION<Type,MType> const&,int);

template <typename Type,typename MType> Type Trace(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> Type Determinant(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> Type Cofactor(MATRIXEXPRESSION<Type,MType> const&,int,int);
template <typename Type,typename MType> Type Spur(MATRIXEXPRESSION<Type,MType> const&);

// *******************************************************************
// ************************* MATRIX INVERSIONS ***********************
// *******************************************************************

template <typename Type,typename MType> MATRIX<Type,0,0> Invert(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> Inverse(MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIX<Type,0,0> DInverse(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> GJInverse(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> LInverse(MATRIXEXPRESSION<Type,MType> const&,bool message=true); // the bool is to whether the warning message is printed
template <typename Type,typename MType> MATRIX<Type,0,0> LTInverse(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> UTInverse(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> LUInverse(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> MPInverse(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> QRInverse(MATRIXEXPRESSION<Type,MType> const&);

// *******************************************************************
// ***************** SPECIAL MATRIX INVERSIONS ***********************
// *******************************************************************

template <typename Type,typename MType> MATRIX<Type,0,0> TridiagonalInverse(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type,typename MType> MATRIX<Type,0,0> VandermondeMatrixInverse(MATRIXEXPRESSION<Type,MType> const&);
template <typename Type> MATRIX<Type,0,0> VandermondeMatrixInverse(std::vector<Type> const&);

// *******************************************************************
// ***************** MATRIX DECOMPOSITIONS ***************************
// *******************************************************************

template <typename MType> std::vector<MATRIX<double,0,0> > LUDecomposition(MATRIXEXPRESSION<double,MType> const&);
template <typename MType> std::vector<MATRIX<std::complex<double>,0,0> > LUDecomposition(MATRIXEXPRESSION<std::complex<double>,MType> const&);

template <typename MType> std::vector<MATRIX<double,0,0> > QRDecomposition(MATRIXEXPRESSION<double,MType> const&);
template <typename MType> std::vector<MATRIX<std::complex<double>,0,0> > QRDecomposition(MATRIXEXPRESSION<std::complex<double>,MType> const&);

template <typename MType> std::vector< MATRIX<double,0,0> > SVDecomposition(MATRIXEXPRESSION<double,MType> const&);
template <typename MType> std::vector< MATRIX<std::complex<double>,0,0> > SVDecomposition(MATRIXEXPRESSION<std::complex<double>,MType> const&);

template <typename Type,typename MType> MATRIX<Type,0,0> Diagonalize(MATRIXEXPRESSION<Type,MType> const&);

// *******************************************************************
// ***************** EIGENVALUES *************************************
// *******************************************************************

template <typename MType> std::vector<std::complex<double> > QREigenValues(MATRIXEXPRESSION<double,MType> const &ME,bool deflate=false);
template <typename MType> std::vector<std::complex<double> > QREigenValues(MATRIXEXPRESSION<std::complex<double>,MType> const &ME,bool deflate=false);

template <typename MType> std::vector<double> JacobiEigenValues(MATRIXEXPRESSION<double,MType> const &ME); // applies to symmetric matrices only
template <typename MType> std::vector<std::complex<double> > JacobiEigenValues(MATRIXEXPRESSION<std::complex<double>,MType> const &ME); // applies to Hermitian matrices only

template <typename Type,typename MType> Type ConditionNumber(MATRIXEXPRESSION<Type,MType> const &ME);
template <typename Type,typename MType> Type InverseConditionNumber(MATRIXEXPRESSION<Type,MType> const &ME);

// *******************************************************************
// ***************** EIGENVECTORS ************************************
// *******************************************************************

template <typename MType> CVECTOR<double,0> EigenVector(MATRIXEXPRESSION<double,MType> const &ME,double Eigenvalue);
template <typename MType> CVECTOR<std::complex<double>,0> EigenVector(MATRIXEXPRESSION<std::complex<double>,MType> const &ME,std::complex<double> Eigenvalue);

// *******************************************************************
// ***************** OTHER MATRIX AND VECTOR OPERATIONS **************
// *******************************************************************

template <typename Type,typename MType> MATRIX<Type,0,0> pow(MATRIXEXPRESSION<Type,MType> const&,int);
template <typename Type,typename MType> MATRIX<Type,0,0> sqrt(MATRIXEXPRESSION<Type,MType> const&);

template <typename Type,typename MType> MATRIX<Type,0,0> SwapRows(MATRIXEXPRESSION<Type,MType> const&,int,int);
template <typename Type,typename MType> MATRIX<Type,0,0> SwapColumns(MATRIXEXPRESSION<Type,MType> const&,int,int);

template <typename Type,typename MType> MATRIX<Type,0,0> Pivot(MATRIXEXPRESSION<Type,MType> const &A,int a,int b);

template <typename MType> MATRIX<double,0,0> SortDiagonalAscending(MATRIXEXPRESSION<double,MType> const&);
template <typename MType> MATRIX<std::complex<double>,0,0> SortDiagonalAscending(MATRIXEXPRESSION<std::complex<double>,MType> const&);

template <typename MType> MATRIX<double,0,0> SortDiagonalDescending(MATRIXEXPRESSION<double,MType> const&);
template <typename MType> MATRIX<std::complex<double>,0,0> SortDiagonalDescending(MATRIXEXPRESSION<std::complex<double>,MType> const&);

template <typename MType> MATRIX<double,0,0> Balance(MATRIXEXPRESSION<double,MType> const&);
template <typename MType> MATRIX<std::complex<double>,0,0> Balance(MATRIXEXPRESSION<std::complex<double>,MType> const&);

template <typename Type,typename MType,typename CVType> CVECTOR<Type,0> GaussEliminate(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y);
template <typename Type,typename MType,typename CVType> CVECTOR<Type,0> DiagonalSolve(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y,bool withtest=true);
template <typename Type,typename MType,typename CVType> CVECTOR<Type,0> TridiagonalSolve(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y,bool withtest=true);
template <typename Type,typename MType,typename CVType> CVECTOR<Type,0> BanddiagonalSolve(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y,int p,int q,bool withtest=true);

template <std::size_t N> CVECTOR<double,N> TridiagonalSolve(MVECTOR<double,N> &L,MVECTOR<double,N> &D,MVECTOR<double,N> &U,CVECTOR<double,N> &Y);
CVECTOR<double,0> TridiagonalSolve(MVECTOR<double,0> &L,MVECTOR<double,0> &D,MVECTOR<double,0> &U,CVECTOR<double,0> &Y);

template <typename Type,typename MType,typename CVType> CVECTOR<Type,0> BackSubstitution(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y,int width,bool withtest=true);
template <typename Type,typename MType,typename CVType> CVECTOR<Type,0> ForwardSubstitution(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y,int width,bool withtest=true);

template <typename MType> MATRIX<double,0,0> UpperHessenberg(MATRIXEXPRESSION<double,MType> const&); // return a matrix with an upper Hessenberg form
template <typename MType> MATRIX<std::complex<double>,0,0> UpperHessenberg(MATRIXEXPRESSION<std::complex<double>,MType> const&); 

template <typename MType> MATRIX<double,0,0> LowerHessenberg(MATRIXEXPRESSION<double,MType> const&); // return a matrix with an lower Hessenberg form
template <typename MType> MATRIX<std::complex<double>,0,0> LowerHessenberg(MATRIXEXPRESSION<std::complex<double>,MType> const&);

template <typename Type,typename MType> CVECTOR<Type,0> Vectorize(MATRIXEXPRESSION<double,MType> const&);
template <typename Type,typename MType> CVECTOR<Type,0> HalfVectorize(MATRIXEXPRESSION<double,MType> const&);

// *******************************************************************
// **************** LINEAR EQUATION SOLVERS **************************
// *******************************************************************

// For when the template type is a double
template <std::size_t N> CVECTOR<double,N> LUSolve(MATRIX<double,N,N> &M,CVECTOR<double,N> &Y);
CVECTOR<double,0> LUSolve(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y);

template <std::size_t N> CVECTOR<double,N> QRSolve(MATRIX<double,N,N> &M,CVECTOR<double,N> &Y);
CVECTOR<double,0> QRSolve(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y);

// *******************************************************************
// *******************************************************************
// *******************************************************************

template <typename Type,typename VType> struct EIGENPAIR
         { Type value;
           VType vector;

           EIGENPAIR(Type VALUE,VType VECTOR) : value(VALUE), vector(VECTOR) {;} 
          };

// *******************************************************************

#endif

