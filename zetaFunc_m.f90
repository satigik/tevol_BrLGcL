module zetaFunc_m
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
  ! !   This is a re-written version of the F77 version, originally written   ! !
  ! !    by "R. L. Mace". Runs slightly faster than the original version.     ! ! 
  ! !                                                                         ! !
  ! !     Evaluates the plasma dispersion function (Fried and Conte           ! !
  ! !     function) of complex argument with a relative error of 1.e-6.       ! !
  ! !                                                                         ! !
  ! !     Algorithm: based closely on that described in Piero Barberio-       ! !
  ! !                Corsetti 'Calculation of the Plasma Dispersion           ! !  
  ! !                Function'.                                               ! ! 
  ! !                                                                         ! ! 
  ! !     Precision: Double                                                   ! ! 
  ! !                                                                         ! !
  ! !     Author: R. L. Mace, Plasma Physics Research Institute               ! !
  ! !               University of Natal, Durban                               ! !
  ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
contains
  complex(wp) function zfn(z)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    complex(wp) :: z
    !! Constants
    real(wp), parameter :: tol=1.e-14_wp,zero=0.e0_wp,half=0.5e0_wp,one=1.e0_wp
    real(wp), parameter :: dlim=6.0_wp,thrhlf=1.5e0_wp,pid4=0.7853981633974480_wp
    complex(wp) :: czero, chalf
    parameter (czero=(0._wp,0._wp), chalf=(0.5_wp,0._wp))
    complex(wp) :: cone, ctwo
    parameter (cone=(1.e0_wp,0.e0_wp), ctwo=(2.e0_wp,0.e0_wp))
    complex(wp) :: irtpi,i2rtpi
    parameter (irtpi=(0.0_wp,1.772453850905516_wp))
    parameter (i2rtpi=(0.0_wp,3.544907701811032_wp))
    !! Local variables
    real(wp) :: x,y,abx,aby,xymax
    real(wp) :: fn,cn,aslim,yasm
    complex(wp) :: an,anm1,bn,bnm1,errz
    complex(wp) :: anp1,bnp1,aa,bb
    complex(wp) :: z2,zinv,summ,term,pterm

    x = real(z)
    y = imag(z)
    abx = abs(x)
    aby = abs(y)
    if(aby > abx) then
       xymax = aby
    else
       xymax = abx
    end if
    fn = zero
    !!     based on the magnitude of the real and imaginary parts of z, 
    !!     determine which of power series, continued fraction, or 
    !!     asymptotic forms to use
    if(aby > one) then
       !! **********************************
       !! employ the continued fraction form
       !! **********************************
       z2 = half-z*z
       an = z
       anm1 = czero
       bn = z2
       bnm1 = cone
       xymax = one/xymax
       !!       compute the continued fraction
       zfn = an/bn
       errz = zfn-anm1/bnm1
       do while(abs(real(errz)) > tol*abs(real(zfn)) .or. &
            abs(imag(errz)) > tol*abs(imag(zfn)))
          fn = fn+one
          cn = xymax/fn
          aa = -fn*(fn-half)*cn
          bb = (z2+fn+fn)*cn
          anp1 = bb*an+aa*anm1
          bnp1 = bb*an+aa*bnm1
          anm1 = an*cn
          an = anp1
          bnm1 = bn*cn
          bn = bnp1
          zfn = an/bn
          errz = zfn-anm1/bnm1
       end do
       
       !!  add the contribution from the pole if Im(z) .le. 0

       if(y < zero) zfn = zfn+i2rtpi*exp(-z*z)
    else if(abx > dlim) then
       !!  ****************************
       !!  use the asmyptotic expansion
       !!  ****************************
       zinv = cone/z
       z2 = chalf*zinv*zinv
       summ = cone
       term = cone
       aslim = x*x+y*y-one
       do while((abs(real(term)) > tol*abs(real(summ)) .or. &
            abs(imag(term)) > tol*abs(imag(summ))) .and. &
            fn <= aslim)
          fn = fn+one
          term = term*(fn+fn-one)*z2
          summ = summ+term
       end do
       zfn = -zinv*summ
       yasm = pid4/abx
       if(y < -yasm) then
          zfn = zfn+i2rtpi*exp(-z*z)
       else if(y <= yasm) then
          zfn = zfn+irtpi*exp(-z*z)
       end if
    else
       !!  *************************
       !!  use the power series form
       !!  ************************* 
       z2 = z*z
       summ = cone
       term = cone
       pterm = cone
       do while(abs(real(term)) > tol*abs(real(summ)) .or. &
            abs(imag(term)) > tol*abs(imag(summ)))
          fn = fn+one
          pterm = pterm*z2/fn
          term = pterm/(fn+fn+one)
          summ = summ+term
       end do
       zfn = (irtpi-ctwo*z*summ)*exp(-z2)
    end if
    return
  end function zfn
end module zetaFunc_m
