#ifndef SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_XERBLA_HPP
#define SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_XERBLA_HPP

#include <iostream>
#include <string>

namespace shape_preserving
{

namespace util { namespace linear_algebra
{

namespace clapack
{

inline int xerbla(const std::string& srname, const int& info,
                  std::ostream& os = std::cout)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    XERBLA  is an error handler for the LAPACK routines.   
    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

    Arguments   
    =========   

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INTEGER   
            The position of the invalid parameter in the parameter list   

            of the calling routine.   

   ===================================================================== 
*/

  os << "** On entry to " << srname
     << ", parameter number " << info
     << " had an illegal value" << std::endl;
  return 0;
} /* xerbla */


} // namespace clapack

}} // namespace util::linear_algebra

} // namespace shape_preserving

#endif // SHAPE_PRESERVING_UTIL_LINEAR_ALGEBRA_CLAPACK_XERBLA_HPP
