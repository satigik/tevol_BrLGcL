module Qsimp_m
contains
  subroutine Qsimp(func,a,b,s)
    implicit none
    ! Version 22/Aug/2015
    real :: a,b,s
    real :: EPS= 1.E-6
    integer :: j
    integer :: jmax= 10
    real :: os,ost,st
    !      EXTERNAL func
    interface 
       real function FUNC(X)
         real, intent(in) :: X
       end function FUNC
    end interface
    ost=-1.E30
    os= -1.E30
    do j=1,JMAX
       call trapzd(func,a,b,st,j)
       s=(4.E0*st-ost)/3.E0
       if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or. &
               (s.eq.0..and.os.eq.0.)) return
       end if
       os=s
       ost=st
    end do
    open(98,file="WARNING_qsimp")
    write(98,*) " Too many steps in qsimp! "
    close(98)
    !      STOP
  end subroutine Qsimp


  subroutine Trapzd(func,a,b,s,n)
    implicit none
    ! Version 22/Aug/2015
    integer :: n
    integer :: it,j
    real :: a,b,s
    real :: del,suma,tnm,x
    !      external func
    interface 
       real function FUNC(X)
         real, intent(in) :: X
       end function FUNC
    end interface
    if (n.eq.1) then
       s=0.5E0*(b-a)*(func(a)+func(b))
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5E0*del
       suma=0.E0
       do j=1,it
          suma=suma+func(x)
          x=x+del
       end do
       s=0.5E0*(s+(b-a)*suma/tnm)
    endif
    return
  end subroutine Trapzd

  subroutine Qsimpb(func,a,b,s)
    implicit none
    ! Version 22/Aug/2015
    ! To be used if "a' is finite, positive, and "b" goes to infinity, or if
    ! "b" is finite, negative, and "a" goes to -infinity.
    real :: a,b,s
    real :: EPS= 1.E-6
    integer :: j
    integer :: jmax= 10
    real :: os,ost,st
    !      EXTERNAL func
    interface 
       real function func(x)
         real, intent(in) :: x
       end function func
    end interface
    ost=-1.E30
    os= -1.E30
    do j=1,JMAX
       call midinf(func,a,b,st,j)
       s=(4.E0*st-ost)/3.E0
       if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or. &
               (s.eq.0.E0.and.os.eq.0.E0)) return
       end if
       os=s
       ost=st
    end do
    open(98,file="WARNING_qsimpb")
    write(98,*) " Too many steps in qsimpb! "
    close(98)
    !      STOP
  end subroutine Qsimpb

  subroutine MIDINF(FUNK,AA,BB,S,N)
    !	Numerical Recipes, Cambridge, 1989, p. 118:
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
    implicit none
    integer :: n
    integer :: it,j
    real :: aa,bb,s,a,b
    real :: del,summ,tnm,x,ddel,func
    !      EXTERNAL func
    interface 
       real function FUNK(X)
         real, intent(in) :: X
       end function FUNK
    end interface
    FUNC(X)=FUNK(1.E0/X)/X**2   ! This is a statement function which effects
    ! the change of variable. 
    B=1.E0/AA 		! These two statements change the limits of
    ! integration accordingly.
    A=1.E0/BB 
    if (N.eq.1) then 		! From this point on, the routine is exactly
       ! identical to MIDPNT.
       S=(B-A)*FUNC(0.5E0*(A+B)) 
       IT=1 
    else 
       TNM=IT 
       DEL=(B-A)/(3.E0*TNM) 
       DDEL=DEL+DEL 
       X=A+0.5E0*DEL 
       SUMM=0.E0 
       ! do 11 J=1,IT
       do J=1,IT   
          SUMM=SUMM+FUNC(X) 
          X=X+DDEL 
          SUMM=SUMM+FUNC(X) 
          X=X+DEL
       end do
       !11     continue 
       S=(S+(B-A)*SUMM/TNM)/3.E0 
       IT=3*IT 
    end if
    return 
  end subroutine MIDINF
end module Qsimp_m
