module DQsimp_m
contains
  subroutine DQsimp(func,a,b,s)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    ! Version 22/Aug/2015
    real(wp) :: a,b,s
    real(wp) :: EPS= 1.E-8_wp
    integer :: j
    integer :: jmax= 25
    real(wp) :: os,ost,st
    !      EXTERNAL func
    interface 
       real(wp) function func(x)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: x
       end function func
    end interface
    ost=-1.e+30_wp
    os= -1.e+30_wp
    do j=1,jmax
       call trapzd(func,a,b,st,j)
       s=(4._wp*st-ost)/3._wp
       if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or. &
               (s.eq.0._wp.and.os.eq.0._wp)) return
       end if
       os=s
       ost=st
    end do
    open(98,file="WARNING_qsimp")
    write(98,*) " Too many steps in qsimp! "
    close(98)
    !      STOP
  end subroutine DQsimp
  subroutine Trapzd(func,a,b,s,n)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    ! Version 22/Aug/2015
    integer :: n
    integer :: it,j
    real(wp) :: a,b,s
    real(wp) :: del,suma,tnm,x
    !      external func
    interface 
       real(wp) function func(x)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: x
       end function func
    end interface
    if (n.eq.1) then
       s=0.5_wp*(b-a)*(func(a)+func(b))
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5_wp*del
       suma=0._wp
       do j=1,it
          suma=suma+func(x)
          x=x+del
       end do
       s=0.5_wp*(s+(b-a)*suma/tnm)
    endif
    return
  end subroutine Trapzd

  subroutine DQsimpb(func,a,b,s)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    ! Version 22/Aug/2015
    ! To be used if "a' is finite, positive, and "b" goes to infinity, or if
    ! "b" is finite, negative, and "a" goes to -infinity.
    real(wp) :: a,b,s
    real(wp) :: EPS= 1.E-8_wp
    integer :: j
    integer :: jmax= 10
    real(wp) :: os,ost,st
    !      EXTERNAL func
    interface 
       real(wp) function func(x)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: x
       end function func
    end interface
    ost=-1.e+30_wp
    os= -1.e+30_wp
    do j=1,JMAX
       call midinf(func,a,b,st,j)
       s=(4._wp*st-ost)/3._wp
       if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or. &
               (s.eq.0._wp.and.os.eq.0._wp)) return
       end if
       os=s
       ost=st
    end do
    open(98,file="WARNING_qsimpb")
    write(98,*) " Too many steps in qsimpb! "
    close(98)
    !      STOP
  end subroutine DQsimpb
  subroutine midinf(funk,aa,bb,s,n)
    !	numerical recipes, Cambridge, 1989, p. 118:
    !	This routine is an exact replacement for MIDPNT, i.e. return as S 
    !	the Nth stage of refinement of the integral of FUNK from AA to BB,
    !	except that the function is evaluated at evenly spaced points in 1/x
    !	rather than in x. This allows the upper limit BB to be as large and
    !	positive as the computer allows, or the lower limit AA to be as large
    !	and negative, but not both. AA and BB must have the same sign.
    ! Version 22/Aug/2015
    !! use precision_m
    ! use cudafor
    ! use openacc
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: n
    integer :: it,j
    real(wp) :: aa,bb,s,a,b
    real(wp) :: del,summ,tnm,x,ddel,func
    !      EXTERNAL func
    interface 
       real(wp) function funk(x)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: x
       end function funk
    end interface
    func(x)=funk(1._wp/x)/x**2   ! this is a statement function which effects
    ! the change of variable. 
    b=1._wp/aa 		! these two statements change the limits of
    ! integration accordingly.
    a=1._wp/bb 
    if (n.eq.1) then 		! from this point on, the routine is exactly
       ! identical to MIDPNT.
       s=(b-a)*func(0.5_wp*(a+b)) 
       it=1 
    else 
       tnm=it 
       del=(b-a)/(3._wp*tnm) 
       ddel=del+del 
       x=a+0.5_wp*del 
       summ=0._wp 
       do 11 j=1,it
       !do j=1,it   
          summ=summ+func(x) 
          x=x+ddel 
          summ=summ+func(x) 
          x=x+del
       !end do
11        continue 
       s=(s+(b-a)*summ/tnm)/3._wp 
       it=3*it 
    end if
    return 
  end subroutine midinf
  
  subroutine DQsimpc(func,a,b,s)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    ! Version 22/Aug/2015
    ! To be used if "a' is finite, positive, and "b" goes to infinity, or if
    ! "b" is finite, negative, and "a" goes to -infinity.
    real(wp) :: a,b,s
    real(wp) :: EPS= 1.E-8_wp
    integer :: j
    integer :: jmax= 20
    real(wp) :: os,ost,st
    !      EXTERNAL func
    interface 
       real(wp) function func(x)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: x
       end function func
    end interface
    ost=-1.e+30_wp
    os= -1.e+30_wp
    do j=1,JMAX
       call midsql(func,a,b,st,j)
       s=(4._wp*st-ost)/3._wp
       if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or. &
               (s.eq.0._wp.and.os.eq.0._wp)) return
       end if
       os=s
       ost=st
    end do
    open(98,file="WARNING_qsimpc")
    write(98,*) " Too many steps in qsimpc! "
    close(98)
    !      STOP
  end subroutine DQsimpc
  subroutine midsql(funk,aa,bb,s,n)
    !	Numerical Recipes, Cambridge, 1989, p. 118:
    !	This routine is an exact replacement for MIDPNT,
    ! except that it allows for an inverse square-root
    ! singularity in the integrand at the lower limit aa .
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: n
    integer :: it,j
    real(wp) :: aa,bb,s,a,b
    real(wp) :: del,summ,tnm
    real(wp) :: x,ddel,func
    
    interface 
       real(wp) function funk(x)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: x
       end function funk
    end interface
    func(x)=2._wp*x*funk(aa+x**2)
    b=sqrt(bb-aa)
    a=0._wp 
    if (n.eq.1) then 		! from this point on, the routine is exactly
       ! identical to midpnt.
       s=(b-a)*func(0.5_wp*(a+b)) 
       it=1 
    else 
       tnm=it 
       del=(b-a)/(3._wp*tnm) 
       ddel=del+del 
       x=a+0.5_wp*del 
       summ=0._wp 
       do 11 j=1,it
          !do j=1,it   
          summ=summ+func(x) 
          x=x+ddel 
          summ=summ+func(x) 
          x=x+del
          !end do
11        continue 
          s=(s+(b-a)*summ/tnm)/3._wp 
          it=3*it 
    end if
    return 
  end subroutine midsql
end module DQsimp_m
