module special_functions
contains
!*****************************************************************************80
!
!! STVH0 computes the Struve function H0(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) SH0, the value of H0(x).
!
function stvh0(x)  result(sh0) 

  implicit none
  real ( kind = 8 ) a0
  real ( kind = 8 ) by0
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km
  real ( kind = 8 ) p0
  real ( kind = 8 ) pi
  real ( kind = 8 ) q0
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) sh0
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  real ( kind = 8 ) ta0
  real ( kind = 8 ) x

  pi = 3.141592653589793D+00
  s = 1.0D+00
  r = 1.0D+00

  if ( x <= 20.0D+00 ) then
    a0 = 2.0D+00 * x / pi
    do k = 1, 60
      r = - r * x / ( 2.0D+00 * k + 1.0D+00 ) * x &
        / ( 2.0D+00 * k + 1.0D+00 )
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
        exit
      end if
    end do

    sh0 = a0 * s

  else

    if ( x < 50.0D+00 ) then
      km = int ( 0.5D+00 * ( x + 1.0D+00 ) )
    else
      km = 25
    end if

    do k = 1, km
      r = - r * ( ( 2.0D+00 * k - 1.0D+00 ) / x ) ** 2
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
        exit
      end if
    end do

    t = 4.0D+00 / x
    t2 = t * t

    p0 = (((( &
      - 0.37043D-05     * t2 &
      + 0.173565D-04 )  * t2 &
      - 0.487613D-04 )  * t2 &
      + 0.17343D-03 )   * t2 &
      - 0.1753062D-02 ) * t2 &
      + 0.3989422793D+00

    q0 = t * ((((( &
        0.32312D-05     * t2 &
      - 0.142078D-04 )  * t2 &
      + 0.342468D-04 )  * t2 &
      - 0.869791D-04 )  * t2 &
      + 0.4564324D-03 ) * t2 &
      - 0.0124669441D+00 )

    ta0 = x - 0.25D+00 * pi
    by0 = 2.0D+00 / sqrt ( x ) &
      * ( p0 * sin ( ta0 ) + q0 * cos ( ta0 ) )
    sh0 = 2.0D+00 / ( pi * x ) * s + by0

  end if

  return
end



!*****************************************************************************80
!
!! STVH1 computes the Struve function H1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) SH1, the value of H1(x).
!
function stvh1(x) result(sh1)

  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) by1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km
  real ( kind = 8 ) p1
  real ( kind = 8 ) pi
  real ( kind = 8 ) q1
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) sh1
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  real ( kind = 8 ) ta1
  real ( kind = 8 ) x

  pi = 3.141592653589793D+00
  r = 1.0D+00

  if ( x <= 20.0D+00 ) then

    s = 0.0D+00
    a0 = - 2.0D+00 / pi
    do k = 1, 60
      r = - r * x * x / ( 4.0D+00 * k * k - 1.0D+00 )
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
        exit
      end if
    end do

    sh1 = a0 * s

  else

    s = 1.0D+00

    if ( x <= 50.0D+00 ) then
      km = int ( 0.5D+00 * x )
    else
      km = 25
    end if

    do k = 1, km
      r = - r * ( 4.0D+00 * k * k - 1.0D+00 ) / ( x * x )
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
        exit
      end if
    end do

    t = 4.0D+00 / x
    t2 = t * t

    p1 = (((( &
        0.42414D-05      * t2 &
      - 0.20092d-04 )    * t2 &
      + 0.580759D-04 )   * t2 &
      - 0.223203D-03 )   * t2 &
      + 0.29218256D-02 ) * t2 &
      + 0.3989422819D+00

    q1 = t * ((((( &
      - 0.36594D-05     * t2 &
      + 0.1622D-04 )    * t2 &
      - 0.398708D-04 )  * t2 &
      + 0.1064741D-03 ) * t2 &
      - 0.63904D-03 )   * t2 &
      + 0.0374008364D+00 )

    ta1 = x - 0.75D+00 * pi
    by1 = 2.0D+00 / sqrt ( x ) * ( p1 * sin ( ta1 ) + q1 * cos ( ta1 ) )
    sh1 = 2.0D+00 / pi * ( 1.0D+00 + s / ( x * x ) ) + by1

  end if

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Compute W=-Ji1=1-xJ0(x)+J1(x)-phi(x), bessel-integral function, first kind, of order 1
!
function W(x) result(y)
real ( kind = 8 ) x
real ( kind = 8 ) y
y =  -x*bessel_j0(x) + bessel_j1(x) - phi(x)+1.d0
end function 


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Compute phi
!
function phi(x) result(y)
real ( kind = 8 ) pi
real ( kind = 8 ) x
real ( kind = 8 ) y

pi = 3.141592653589793D+00
y = pi/2*x*(bessel_j1(x)*stvh0(x)-bessel_j0(x)*stvh1(x))
end function 



end module







