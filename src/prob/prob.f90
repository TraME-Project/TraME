!
!! Modified by Keith O'Hara
!
!     27 May 2016
!
!! Removed 'timestamp' subroutine
!

subroutine beta_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! BETA_CDF evaluates the Beta CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta_inc
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    cdf = 0.0D+00
  else if ( x <= 1.0D+00 ) then
    cdf = beta_inc ( a, b, x )
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine beta_cdf_inv ( cdf, p, q, x )

!*****************************************************************************80
!
!! BETA_CDF_INV computes the inverse of the incomplete Beta function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    GW Cran, KJ Martin, GE Thomas,
!    Remark AS R19 and Algorithm AS 109:
!    A Remark on Algorithms AS 63: The Incomplete Beta Integral
!    and AS 64: Inverse of the Incomplete Beta Integeral,
!    Applied Statistics,
!    Volume 26, Number 1, 1977, pages 111-114.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the Beta CDF.
!    0 <= CDF <= 1.
!
!    Input, real ( kind = 8 ) P, Q, the parameters of the incomplete
!    Beta function.
!
!    Output, real ( kind = 8 ) X, the argument of the incomplete
!    Beta function which produces the value CDF.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) SAE, the most negative decimal exponent
!    which does not cause an underflow.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) acu
  real ( kind = 8 ) adj
  real ( kind = 8 ) beta_inc
  real ( kind = 8 ) beta_log
  real ( kind = 8 ) cdf
  real ( kind = 8 ) fpu
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) iex
  logical indx
  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ) prev
  real ( kind = 8 ) q
  real ( kind = 8 ) qq
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ), parameter :: sae = -37.0D+00
  real ( kind = 8 ) sq
  real ( kind = 8 ) t
  real ( kind = 8 ) tx
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xin
  real ( kind = 8 ) y
  real ( kind = 8 ) yprev

  fpu = 10.0D+00 ** sae
  beta_log = lgamma ( p ) + lgamma ( q ) - lgamma ( p + q )
!
!  Test for admissibility of parameters.
!
  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  P <= 0.0'
    return
  end if

  if ( q <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  Q <= 0.0'
    return
  end if

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0.0 or 1.0 < CDF.'
    return
  end if
!
!  Return immediately if the answer is easy to determine.  
!
  if ( cdf == 0.0D+00  ) then
    x = 0.0D+00
    return
  else if ( cdf == 1.0D+00 ) then
    x = 1.0D+00
    return
  end if
!
!  Change tail if necessary.
!
  if ( 0.5D+00 < cdf ) then
    a = 1.0D+00 - cdf
    pp = q
    qq = p
    indx = .true.
  else
    a = cdf
    pp = p
    qq = q
    indx = .false.
  end if
!
!  Calculate the initial approximation.
!
  r = sqrt ( - log ( a * a ) )

  y = r - ( 2.30753D+00 + 0.27061D+00 * r ) &
    / ( 1.0D+00 + ( 0.99229D+00 + 0.04481D+00 * r ) * r )

  if ( 1.0D+00 < pp .and. 1.0D+00 < qq ) then

    r = ( y * y - 3.0D+00 ) / 6.0D+00
    s = 1.0D+00 / ( pp + pp - 1.0D+00 )
    t = 1.0D+00 / ( qq + qq - 1.0D+00 )
    h = 2.0D+00 / ( s + t )
    w = y * sqrt ( h + r ) / h - ( t - s ) &
    * ( r + 5.0D+00 / 6.0D+00 - 2.0D+00 / ( 3.0D+00 * h ) )
    x = pp / ( pp + qq * exp ( w + w ) )

  else

    r = qq + qq
    t = 1.0D+00 / ( 9.0D+00 * qq )
    t = r * ( 1.0D+00 - t + y * sqrt ( t ) ) ** 3

    if ( t <= 0.0D+00 ) then
      x = 1.0D+00 - exp ( ( log ( ( 1.0D+00 - a ) * qq ) &
        + beta_log ) / qq )
    else

      t = ( 4.0D+00 * pp + r - 2.0D+00 ) / t

      if ( t <= 1.0D+00 ) then
        x = exp ( ( log ( a * pp ) + beta_log ) / pp )
      else
        x = 1.0D+00 - 2.0D+00 / ( t + 1.0D+00 )
      end if

    end if

  end if
!
!  Solve for X by a modified Newton-Raphson method.
!
  r = 1.0D+00 - pp
  t = 1.0D+00 - qq
  yprev = 0.0D+00
  sq = 1.0D+00
  prev = 1.0D+00

  if ( x < 0.0001D+00 ) then
    x = 0.0001D+00
  end if

  if ( 0.9999D+00 < x ) then
    x = 0.9999D+00
  end if

  iex = max ( - 5.0D+00 / pp ** 2 - 1.0D+00 / a ** 0.2D+00 - 13.0D+00, sae )

  acu = 10.0D+00 ** iex

  do

    y = beta_inc ( pp, qq, x )

    xin = x
    y = ( y - a ) * exp ( beta_log + r * log ( xin ) &
      + t * log ( 1.0D+00 - xin ) )

    if ( y * yprev <= 0.0D+00 ) then
      prev = max ( sq, fpu )
    end if

    g = 1.0D+00

    do

      do

        adj = g * y
        sq = adj * adj

        if ( sq < prev ) then

          tx = x - adj

          if ( 0.0D+00 <= tx .and. tx <= 1.0D+00 ) then
            exit
          end if

        end if

        g = g / 3.0D+00

      end do

      if ( prev <= acu ) then
        if ( indx ) then
          x = 1.0D+00 - x
        end if
        return
      end if

      if ( y * y <= acu ) then
        if ( indx ) then
          x = 1.0D+00 - x
        end if
        return
      end if

      if ( tx /= 0.0D+00 .and. tx /= 1.0D+00 ) then
        exit
      end if

      g = g / 3.0D+00

    end do

    if ( tx == x ) then
      exit
    end if

    x = tx
    yprev = y

  end do

  if ( indx ) then
    x = 1.0D+00 - x
  end if

  return
end

subroutine beta_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_CDF_VALUES returns some values of the Beta CDF.
!
!  Discussion:
!
!    The incomplete Beta function may be written
!
!      BETA_INC(A,B,X) = integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
!                      / integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
!
!    Thus,
!
!      BETA_INC(A,B,0.0) = 0.0
!      BETA_INC(A,B,1.0) = 1.0
!
!    The incomplete Beta function is also sometimes called the
!    "modified" Beta function, or the "normalized" Beta function
!    or the Beta CDF (cumulative density function).
!
!    In Mathematica, the function can be evaluated by:
!
!      BETA[X,A,B] / BETA[A,B]
!
!    The function can also be evaluated by using the Statistics package:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = BetaDistribution [ a, b ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968,
!    LC: QA351.P38.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 45

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.5D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    30.0D+00, &
    30.0D+00, &
    40.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     1.30625D+00, &
     1.30625D+00, &
     1.30625D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.0D+00, &
     0.5D+00, &
     5.0D+00, &
     5.0D+00, &
    10.0D+00, &
     5.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
    11.7562D+00, &
    11.7562D+00, & 
    11.7562D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6376856085851985D-01, &
    0.2048327646991335D+00, &
    0.1000000000000000D+01, &
    0.0000000000000000D+00, &
    0.5012562893380045D-02, &
    0.5131670194948620D-01, &
    0.2928932188134525D+00, &
    0.5000000000000000D+00, &
    0.2800000000000000D-01, &
    0.1040000000000000D+00, &
    0.2160000000000000D+00, &
    0.3520000000000000D+00, &
    0.5000000000000000D+00, &
    0.6480000000000000D+00, &
    0.7840000000000000D+00, &
    0.8960000000000000D+00, &
    0.9720000000000000D+00, &
    0.4361908850559777D+00, &
    0.1516409096347099D+00, &
    0.8978271484375000D-01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.4598773297575791D+00, &
    0.2146816102371739D+00, &
    0.9507364826957875D+00, &
    0.5000000000000000D+00, &
    0.8979413687105918D+00, &
    0.2241297491808366D+00, &
    0.7586405487192086D+00, &
    0.7001783247477069D+00, &
    0.5131670194948620D-01, &
    0.1055728090000841D+00, &
    0.1633399734659245D+00, &
    0.2254033307585166D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.6723200000000000D+00, &
    0.2160000000000000D+00, &
    0.8370000000000000D-01, &
    0.3078000000000000D-01, &
    0.1093500000000000D-01, &
    0.918884684620518D+00, &
    0.21052977489419D+00, &
    0.1824130512500673D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.10D+00, &
    1.00D+00, &
    0.00D+00, &
    0.01D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    0.50D+00, &
    0.90D+00, &
    0.50D+00, &
    1.00D+00, &
    0.50D+00, &
    0.80D+00, &
    0.60D+00, &
    0.80D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.70D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.225609D+00, &
    0.0335568D+00, &
    0.0295222D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function beta_check ( a, b )

!*****************************************************************************80
!
!! BETA_CHECK checks the parameters of the Beta PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, logical BETA_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical beta_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    beta_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    beta_check = .false.
    return
  end if

  beta_check = .true.

  return
end
function beta_inc ( a, b, x )

!*****************************************************************************80
!
!! BETA_INC returns the value of the incomplete Beta function.
!
!  Discussion:
!
!    This calculation requires an iteration.  In some cases, the iteration
!    may not converge rapidly, or may become inaccurate.
!
!    The formula is:
!
!      BETA_INC(A,B,X)
!
!        =   Integral ( 0 <= T <= X ) T^(A-1) (1-T)^(B-1) dT
!          / Integral ( 0 <= T <= 1 ) T^(A-1) (1-T)^(B-1) dT
!
!        =   Integral ( 0 <= T <= X ) T^(A-1) (1-T)^(B-1) dT
!          / BETA(A,B)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    KL Majumder, GP Bhattacharjee,
!    Algorithm AS63,
!    Applied Statistics,
!    1973, volume 22, number 3.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the function.
!    0.0 < A,
!    0.0 < B.
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!    Normally, 0.0D+00 <= X <= 1.0.
!
!    Output, real ( kind = 8 ) BETA_INC, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta_inc
  real ( kind = 8 ) cx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 1000
  logical indx
  integer ( kind = 4 ) ns
  real ( kind = 8 ) pp
  real ( kind = 8 ) psq
  real ( kind = 8 ) qq
  real ( kind = 8 ) r8_beta
  real ( kind = 8 ) rx
  real ( kind = 8 ) temp
  real ( kind = 8 ) term
  real ( kind = 8 ), parameter :: tol = 1.0D-07
  real ( kind = 8 ) x
  real ( kind = 8 ) xx

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_INC - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    stop 1
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_INC - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    stop 1
  end if

  if ( x <= 0.0D+00 ) then
    beta_inc = 0.0D+00
    return
  else if ( 1.0D+00 <= x ) then
    beta_inc = 1.0D+00
    return
  end if
!
!  Change tail if necessary and determine S.
!
  psq = a + b

  if ( a < ( a + b ) * x ) then
    xx = 1.0D+00 - x
    cx = x
    pp = b
    qq = a
    indx = .true.
  else
    xx = x
    cx = 1.0D+00 - x
    pp = a
    qq = b
    indx = .false.
  end if

  term = 1.0D+00
  i = 1
  beta_inc = 1.0D+00

  ns = int ( qq + cx * ( a + b ) )
!
!  Use Soper's reduction formulas.
!
  rx = xx / cx

  temp = qq - real ( i, kind = 8 )
  if ( ns == 0 ) then
    rx = xx
  end if

  it = 0

  do

    it = it + 1

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BETA_INC - Fatal error!'
      write ( *, '(a)' ) '  Maximum number of iterations exceeded!'
      write ( *, '(a,i8)' ) '  IT_MAX = ', it_max
      stop 1
    end if

    term = term * temp * rx / ( pp + real ( i, kind = 8 ) )
    beta_inc = beta_inc + term
    temp = abs ( term )

    if ( temp <= tol .and. temp <= tol * beta_inc ) then
      exit
    end if

    i = i + 1
    ns = ns - 1

    if ( 0 <= ns ) then
      temp = qq - real ( i, kind = 8 )
      if ( ns == 0 ) then
        rx = xx
      end if
    else
      temp = psq
      psq = psq + 1.0D+00
    end if

  end do
!
!  Finish calculation.
!
  beta_inc = beta_inc * exp ( pp * log ( xx ) &
    + ( qq - 1.0D+00 ) * log ( cx ) ) / ( r8_beta ( a, b ) * pp )

  if ( indx ) then
    beta_inc = 1.0D+00 - beta_inc
  end if

  return
end
subroutine beta_inc_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_INC_VALUES returns some values of the incomplete Beta function.
!
!  Discussion:
!
!    The incomplete Beta function may be written
!
!      BETA_INC(A,B,X) = integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
!                      / integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
!
!    Thus,
!
!      BETA_INC(A,B,0.0) = 0.0
!      BETA_INC(A,B,1.0) = 1.0
!
!    The incomplete Beta function is also sometimes called the
!    "modified" Beta function, or the "normalized" Beta function
!    or the Beta CDF (cumulative density function).
!
!    In Mathematica, the function can be evaluated by:
!
!      BETA[X,A,B] / BETA[A,B]
!
!    The function can also be evaluated by using the Statistics package:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = BetaDistribution [ a, b ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968,
!    LC: QA351.P38.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 45

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.5D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    30.0D+00, &
    30.0D+00, &
    40.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     1.30625D+00, &
     1.30625D+00, &
     1.30625D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.0D+00, &
     0.5D+00, &
     5.0D+00, &
     5.0D+00, &
    10.0D+00, &
     5.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
    11.7562D+00, &
    11.7562D+00, & 
    11.7562D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6376856085851985D-01, &
    0.2048327646991335D+00, &
    0.1000000000000000D+01, &
    0.0000000000000000D+00, &
    0.5012562893380045D-02, &
    0.5131670194948620D-01, &
    0.2928932188134525D+00, &
    0.5000000000000000D+00, &
    0.2800000000000000D-01, &
    0.1040000000000000D+00, &
    0.2160000000000000D+00, &
    0.3520000000000000D+00, &
    0.5000000000000000D+00, &
    0.6480000000000000D+00, &
    0.7840000000000000D+00, &
    0.8960000000000000D+00, &
    0.9720000000000000D+00, &
    0.4361908850559777D+00, &
    0.1516409096347099D+00, &
    0.8978271484375000D-01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.4598773297575791D+00, &
    0.2146816102371739D+00, &
    0.9507364826957875D+00, &
    0.5000000000000000D+00, &
    0.8979413687105918D+00, &
    0.2241297491808366D+00, &
    0.7586405487192086D+00, &
    0.7001783247477069D+00, &
    0.5131670194948620D-01, &
    0.1055728090000841D+00, &
    0.1633399734659245D+00, &
    0.2254033307585166D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.6723200000000000D+00, &
    0.2160000000000000D+00, &
    0.8370000000000000D-01, &
    0.3078000000000000D-01, &
    0.1093500000000000D-01, &
    0.918884684620518D+00, &
    0.21052977489419D+00, &
    0.1824130512500673D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.10D+00, &
    1.00D+00, &
    0.00D+00, &
    0.01D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    0.50D+00, &
    0.90D+00, &
    0.50D+00, &
    1.00D+00, &
    0.50D+00, &
    0.80D+00, &
    0.60D+00, &
    0.80D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.70D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.225609D+00, &
    0.0335568D+00, &
    0.0295222D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end

subroutine beta_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! BETA_PDF evaluates the Beta PDF.
!
!  Discussion:
!
!    The formula for the PDF is:
!
!      PDF(A,B;X) = X^(A-1) * (1-X)^(B-1) / BETA(A,B).
!
!    A = B = 1 yields the Uniform distribution on [0,1].
!    A = B = 1/2 yields the Arcsin distribution.
!        B = 1 yields the power function distribution.
!    A = B -> Infinity tends to the Normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0 <= X <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) r8_beta
  real ( kind = 8 ) x

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    pdf = 0.0D+00
  else
    pdf = x ** ( a - 1.0D+00 ) * ( 1.0D+00 - x ) ** ( b - 1.0D+00 ) &
      / r8_beta ( a, b )
  end if

  return
end

subroutine beta_values ( n_data, x, y, fxy )

!*****************************************************************************80
!
!! BETA_VALUES returns some values of the Beta function.
!
!  Discussion:
!
!    Beta(X,Y) = ( Gamma(X) * Gamma(Y) ) / Gamma(X+Y)
!
!    Both X and Y must be greater than 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      Beta[X,Y]
!
!  Properties:
!
!    Beta(X,Y) = Beta(Y,X).
!    Beta(X,Y) = integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
!    Beta(X,Y) = Gamma(X) * Gamma(Y) / Gamma(X+Y)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, Y, the arguments of the function.
!
!    Output, real ( kind = 8 ) FXY, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    0.5000000000000000D+01, &
    0.2500000000000000D+01, &
    0.1666666666666667D+01, &
    0.1250000000000000D+01, &
    0.5000000000000000D+01, &
    0.2500000000000000D+01, &
    0.1000000000000000D+01, &
    0.1666666666666667D+00, &
    0.3333333333333333D-01, &
    0.7142857142857143D-02, &
    0.1587301587301587D-02, &
    0.2380952380952381D-01, &
    0.5952380952380952D-02, &
    0.1984126984126984D-02, &
    0.7936507936507937D-03, &
    0.3607503607503608D-03, &
    0.8325008325008325D-04 /)
  real ( kind = 8 ) fxy
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    6.0D+00, &
    7.0D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( n_max ) :: y_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    0.2D+00, &
    0.4D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    5.0D+00, &
    6.0D+00, &
    7.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    y = 0.0D+00
    fxy = 0.0D+00
  else
    x = x_vec(n_data)
    y = y_vec(n_data)
    fxy = b_vec(n_data)
  end if

  return
end

function r8_beta ( a, b )

!*****************************************************************************80
!
!! R8_BETA returns the value of the Beta function.
!
!  Discussion:
!
!    The Beta function is defined as
!
!      BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
!                = Integral ( 0 <= T <= 1 ) T^(A-1) (1-T)^(B-1) dT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the function.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real ( kind = 8 ) R8_BETA, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r8_beta

  if ( a <= 0.0D+00 .or. b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BETA - Fatal error!'
    write ( *, '(a)' ) '  Both A and B must be greater than 0.'
    stop 1
  end if

  r8_beta = exp ( lgamma ( a ) + lgamma ( b ) - lgamma ( a + b ) )

  return
end
