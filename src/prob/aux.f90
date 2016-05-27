!
!! Auxiliary functions for TraME
!
!     27 May 2016
!

subroutine beta_cdf_inv_int ( x, p, q, res )

!*****************************************************************************80
!
!! BETA_CDF_INV_INT integrates the beta dist. quantile function using DQAG.
!
!  Discussion:
!
!    Uses default options for DQAG
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    Keith O'Hara
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) res

  integer ( kind = 4 ), parameter :: limit = 500
  integer ( kind = 4 ), parameter :: lenw = 4 * limit

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: beta_cdf_inv_fn
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ), parameter :: key = 6
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) work(lenw)
  
  real ( kind = 8 ) f_args(2)
  f_args(1) = p
  f_args(2) = q

  b = x

  call dqag ( beta_cdf_inv_fn, f_args, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
    limit, lenw, last, iwork, work )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test beta_cdf_inv_int'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  res = result
  
  return
end


function beta_cdf_inv_fn ( cdf, f_args )

!*****************************************************************************80
!
!! BETA_CDF_INV_FN computes the inverse of the incomplete Beta function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    Keith O'Hara
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the Beta CDF.
!    0 <= CDF <= 1.
!
!    Input, real ( kind = 8 ) f_args, the parameters of the incomplete
!    Beta function.
!
!    Output, real ( kind = 8 ) beta_cdf_inv_fn, the argument of the incomplete
!    Beta function which produces the value CDF.
!

  real ( kind = 8 ) cdf
  real ( kind = 8 ) f_args(2)
  
  real ( kind = 8 ) beta_cdf_inv_fn
  
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  
  p = f_args(1)
  q = f_args(2)
  
  call beta_cdf_inv(cdf,p,q,beta_cdf_inv_fn)
  
  return
end