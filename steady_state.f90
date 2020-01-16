  
module common_params
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp) :: Qxi=1.e-4_wp, Qxf=0.6_wp
  real(wp) :: Qzi=1.e-4_wp, Qzf=0.6_wp
  real(wp) :: Qi= 1.e-4_wp, Qf= 6.0_wp ! For VQQ
  real(wp) :: Ulim= 16.0_wp
  real(wp) :: RTeTi= 1._wp, G= 5.e-3_wp
  real(wp) :: Ve2C2= 4.e-3_wp
  real(wp) :: RMeMi, RMiMe
  real(wp) :: Geff, AA
  real(wp) :: kappa= 2.5_wp, NkNe= 0.1_wp
  character(LEN=3) :: BremGcollL= "Yes"
  integer :: nqx=64, nqz=64
  integer :: nux=64, nuz=128
  integer :: nqcd=1024
  integer :: nrz=1
end module common_params

module common_arrays
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp), dimension(:), allocatable :: VQ,VQQ
  real(wp), dimension(:), allocatable :: VQx
  real(wp), dimension(:), allocatable :: VQz
  real(wp), dimension(:), allocatable :: VUx,VUz
  real(wp), dimension(:,:), allocatable :: Fe,Fe0,Fek
  real(wp), dimension(:), allocatable :: Fe1D,Fe01D,Fek1D
  real(wp), dimension(:,:), allocatable :: IL,IL0
  real(wp), dimension(:,:), allocatable :: ILk,ILMaxw
  real(wp), dimension(:), allocatable :: IL1D,IL01D
  real(wp), dimension(:), allocatable :: ILk1D,ILMaxw1D
  real(wp), dimension(:,:), allocatable :: GcollLp
  real(wp), dimension(:,:), allocatable :: GcollLm
  ! real(wp), dimension(:,:), allocatable :: GcollL 
  real(wp), dimension(:), allocatable:: GcollL1D
  real(wp), dimension(:), allocatable:: GcollLp1D,GcollLm1D
  ! real(wp), dimension(:,:), allocatable :: BremL
  real(wp), dimension(:,:), allocatable :: BremLp,BremLm
  real(wp), dimension(:), allocatable :: BremL1D
  ! real(wp), dimension(:,:), allocatable :: BrGcL
  ! real(wp), dimension(:), allocatable :: BrGcL1D
  real(wp), dimension(:,:), allocatable :: Vint
  real(wp), dimension(:), allocatable :: Vint1D
  ! real(wp), dimension(:,:), allocatable :: GqlLp,GqlLm
  ! real(wp), dimension(:,:), allocatable :: GqlL1D
  ! real(wp), dimension(:,:), allocatable:: GcollSp,GqlSp
  ! real(wp), dimension(:,:), allocatable:: GcollSm,GqlSm
  ! real(wp), dimension(:), allocatable:: GcollS1D,GqlS1D,BremS1D
  ! real(wp), dimension(:), allocatable:: GcollSp1D,GqlSp1D
  ! real(wp), dimension(:), allocatable:: GcollSm1D,GqlSm1D
  ! real(wp), dimension(:,:), allocatable :: BremSp,BremSm,BremS
  real(wp), dimension(:), allocatable :: VRNeNs,VRTeTs,VRTiTs
  ! Auxiliary for collisional damping:
  ! real(wp), dimension(:), allocatable :: Aux1_Gcoll
  ! integer, dimension(:), allocatable :: Aux2_Gcoll
  ! Auxiliary for electrostatic bremsstrahlung:
  ! real(wp), dimension(:), allocatable :: Aux1_Bremss
  ! integer, dimension(:), allocatable :: Aux2_Bremss
end module common_arrays

module math_constants 
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp), parameter :: Pi = 4._wp*atan(1._wp)
  real(wp), parameter :: Sqtwo = sqrt(2._wp)
  ! real(wp) :: Pi= 3.1415926535897932384626433832795_wp
  ! real(wp) :: Sqtwo= 1.41421356237309504880168872420969_wp
  real(wp), parameter :: Infinity = 1.0e+30_wp
  real(wp), parameter :: Degree = 0.01745329251994329576923690768488_wp
  real(wp), parameter :: EpsMin = 1.0e-16_wp
  real(wp), parameter :: Xacc = 1.0e-6_wp
  real(wp), parameter :: Qmin = 1.0e-4_wp
  real(wp), parameter :: Dqaux = 1.0e-2_wp
  complex(wp), parameter :: Zi = (0._wp,1._wp)
  complex(wp), parameter :: Zzero = (0._wp,0._wp)
end module math_constants

module phys_constants
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp), parameter :: Me = 9.10938356e-28_wp       ! electron mass (g)
  real(wp), parameter :: Mi = 1.6726219e-24_wp        ! proton mass (g)
  real(wp), parameter :: MeC2 = 510.998946e+3_wp      ! electron mass (eV)
  real(wp), parameter :: MpC2 = 938.272081e+6_wp      ! proton mass (eV)
  real(wp), parameter :: C_SI = 2.997925e+8_wp        ! speed of light (m/s)
  real(wp), parameter :: C_cgs = 2.997925e+10_wp      ! speed of light (cm/s)  
end module phys_constants

module sub_prog
contains

  subroutine allocate_arrays
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    implicit none
    allocate (VQ(nrz))
    allocate (VQx(nqx))
    allocate (VQz(nqz))
    allocate (VQQ(nqcd))
    allocate (VUx(nux),VUz(nuz))
    allocate (Fe(nux,nuz),Fe0(nux,nuz))
    allocate (Fek(nux,nuz))
    allocate (Fe1D(nuz),Fe01D(nuz))
    allocate (Fek1D(nuz))
    allocate (IL0(nqx,nqz),IL(nqx,nqz))
    allocate (ILk(nqx,nqz),ILMaxw(nqx,nqz))
    allocate (IL01D(nqcd),IL1D(nqcd))
    allocate (ILk1D(nqcd),ILMaxw1D(nqcd))
    allocate (GcollLp(nqx,nqz),GcollLm(nqx,nqz))
    ! allocate (GcollL(nqx,nqz))
    allocate (GcollL1D(nqcd))
    allocate (GcollLp1D(nqcd),GcollLm1D(nqcd))
    ! allocate (BremL(nqx,nqz))
    allocate (BremL1D(nqcd))
    allocate (BremLp(nqx,nqz))
    allocate (BremLm(nqx,nqz))
    ! allocate (BrGcL(nqx,nqz))
    ! allocate (BrGcL1D(nqcd))
    allocate (Vint(nqx,nqz))
    allocate (Vint1D(nqcd))
    ! allocate (GqlLp(nqx,nqz),GqlLm(nqx,nqz))
    ! allocate (GqlL1D(nqcd))
    allocate (VRNeNs(nrz),VRTeTs(nrz),VRTiTs(nrz))
    ! Auxiliary for collisional damping:
    ! allocate (Aux1_Gcoll(4))
    ! allocate (Aux2_Gcoll(2))
    ! Auxiliary for electrostatic bremsstrahlung:
    ! allocate(Aux1_Bremss(5))
    ! allocate(Aux2_bremss(2))
  end subroutine allocate_arrays

  subroutine Definitions
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    use math_constants
    use phys_constants
    implicit none
    real(wp) :: RMpMe,RMeMp
    integer :: m

    m=1
    RMpMe= MpC2/MeC2
    RMeMp= MeC2/MpC2
    RMiMe= RMpMe
    RMeMi= 1._wp/RMiMe
    AA= sqrt(1._wp+3._wp/RTeTi)/sqrt(RMiMe)/sqrt(2._wp)
    !Geff= G/(2._wp*sqrt(2._wp)*(4._wp*Pi)**2)
    VRNeNs(m)= 1._wp
    VRTeTs(m)= 1._wp
    VRTiTs(m)= 1._wp/RTeTi
  end subroutine Definitions

  subroutine init_vec
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    implicit none
    real(wp) :: Ux,Uz
    real(wp) :: Dqx,Dqz,Dq
    real(wp) :: Dux,Duz
    ! real(wp) :: Dux2,Duz2
    ! real(wp), parameter :: Ucrit=1.6_wp
    integer :: i,nu2

    !! !Wave-vector space
    Dqx= (Qxf-Qxi)/(nqx-1)
    do i= 1,nqx
       VQx(i)= Qxi+(i-1)*Dqx
    end do
    VQx(nqx)= Qxf
    
    Dqz= (Qzf-Qzi)/(nqz-1)
    do i= 1,nqz
       VQz(i)= Qzi+(i-1)*Dqz
    end do
    VQz(nqz)= Qzf

    Dq= (Qf-Qi)/(nqcd-1)
    do i= 1,nqcd
       VQQ(i)= Qi+(i-1)*Dq
    end do
    VQQ(nqcd)= Qf

    !! !Velocity space
    Dux= (Ulim-0._wp)/(nux-1)
    do i= 1,nux
       VUx(i)= 0._wp+(i-1)*Dux
    end do
    
    Duz= (Ulim-(-Ulim))/(nuz-1)    
    do i= 1,nuz
       VUz(i)= -Ulim+(i-1)*Duz
    end do

    !! !Symmetrization of the array for Uz:
    if(mod(nuz,2) .ne. 0) then
       nu2= (nuz-1)/2+1
       do i= 1,nu2
          VUz(nuz+1-i)= -VUz(i)
       end do
       VUz(nu2)= 0._wp
    end if
  end subroutine init_vec

  subroutine fe_init
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use common_arrays
    use common_params
    use math_constants
    implicit none
    real(wp) :: Ux,Uz
    real(wp) :: U,U2
    integer :: i,k,m
    
    m=1
    
    do k= 1,nuz
       Uz= VUz(k)
       do i= 1,nux
          Ux= VUx(i)
          U2= Ux**2+Uz**2
          U= sqrt(U2)
          !! !Maxwellian
          Fe0(i,k)= 1._wp/(Pi*VRTeTs(m))*exp(-U2/VRTeTs(m))
          ! !! Core-halo EVDF! !!
          ! Fe(i,k)= (1._wp-NkNe)/Pi/VRTeTs(m)*exp(-U2/VRTeTs(m)) &
          !      + NkNe/Pi/VRTeTs(m)*kappa/(kappa-1._wp) &
          !      * (1._wp+U2/(kappa-1._wp)/VRTeTs(m))**(-(kappa+1._wp))
          ! !! !Kappa EVDF! !!
          ! Fek(i,k)= 1._wp/Pi/VRTeTs(m)*kappa/(kappa-1._wp) &
          !      * (1._wp+U2/(kappa-1._wp)/VRTeTs(m))**(-(kappa+1._wp))
       end do
    end do
    
    call output("Fe ")
    return
  end subroutine fe_init

  subroutine wave_init
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use common_arrays
    use common_params
    use math_constants
    use integralBrGc_m
    implicit none
    real(wp) :: Q,Q2
    real(wp) :: Qx,Qz
    real(wp) :: Zlq,Zlq2
    real(wp) :: Aux
    integer :: i,j,m,ires
    m = 1
    Geff= G/(2._wp*sqrt(2._wp)*(4._wp*pi)**2)

    select case(BremGcollL)
    case("No ")
       do i=1,nqcd
          Q= VQQ(i)
          Q2=Q*Q
          Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
          Zlq2= Zlq*Zlq
          !! !Initial spectrum for L waves for Maxwellian EVDF! !!
          IL01D(i)= Geff*VRNeNs(m)*VRTeTs(m)/2._wp/Zlq2
          !! !Initial spectrum for L waves for core-halo EVDF! !!
          ! IL1D(i)= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m))/Q/Q2 &
          !      * ((1._wp-NkNe)*exp(-Zlq2/Q2/VRTeTs(m)) &
          !      + NkNe/sqrt(kappa-1.5_wp)*gamma(kappa)/gamma(kappa-0.5_wp) &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(kappa-1.5_wp))**(-kappa))) &
          !      / (2._wp*sqrt(Pi)*VRNeNS(m)*(Zlq2)/sqrt(VRTeTs(m))/Q/Q2 &
          !      * ((1._wp-NkNe)*exp(-Zlq2/Q2/VRTeTs(m)) &
          !      + NkNe/sqrt((kappa-1.5_wp)**3)*gamma(kappa+1._wp)/gamma(kappa-0.5_wp) &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(kappa-1.5_wp))**(-(kappa+1._wp))))          
          ! !! !Initial spectrum of L waves for kappa EVDF! !!
          ! ILk1D(i)= Geff*VRNeNs(m)*VRTeTs(m)*(kappa-1.5_wp)/kappa/Zlq2/2._wp &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(kappa-1.5_wp))
       end do
       Q= 0._wp
       Q2= 0._wp
       Zlq= 0._wp
       Zlq2= 0._wp
       do i=1,nqx
          Qx=VQx(i)
          do j=1,nqz
             Qz=VQz(j)
             Q2=Qx**2+Qz**2
             Q=sqrt(Q2)
             ! call Locate(VQQ,nqcd,Q,ires)
             ! call Aitp1d2(nqcd,VQQ,IL1D,Q,Aux,ires)
             ! IL(i,j)=Aux
             ! call Locate(VQQ,nqcd,Q,ires)
             ! call Aitp1d2(nqcd,VQQ,ILk1D,Q,Aux,ires)
             ! ILk(i,j)=Aux
             Zlq= ZL(Qx,Qz)
             Zlq2= Zlq**2
             IL0(i,j)=Geff*VRNeNs(m)*VRTeTs(m)/2._wp/Zlq2
          end do
       end do
       call output("IL ")
    case("Yes")
       call integral(nqcd,VQQ,RTeTi,Vint1D)
       call coll_damping
       call bremsstrahlung
       do i=1,nqcd
          Q= VQQ(i)
          Q2=Q*Q
          Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
          Zlq2= Zlq*Zlq
          !! !Only spontt. and ind. emissions with Maxw EVDF! !!
          IL01D(i)= Geff*VRNeNs(m)*VRTeTs(m)/2._wp/Zlq2
          
          ! !! !Core-halo EVDF with spont. and ind. emissions, including EB and CD
          ! IL1D(i)= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m))/Q/Q2 &
          !      * ((1._wp-NkNe)*exp(-Zlq2/Q2/VRTeTs(m)) &
          !      + NkNe/sqrt(kappa-1.5_wp)*gamma(kappa)/gamma(kappa-0.5_wp) &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(kappa-1.5_wp))**(-kappa))+BremL1D(i)) &
          !      / (2._wp*sqrt(Pi)*VRNeNS(m)*(Zlq2)/sqrt(VRTeTs(m))/Q/Q2 &
          !      * ((1._wp-NkNe)*exp(-Zlq2/Q2/VRTeTs(m)) &
          !      + NkNe/sqrt((kappa-1.5_wp)**3)*gamma(kappa+1._wp)/gamma(kappa-0.5_wp) &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(kappa-1.5_wp))**(-(kappa+1._wp)))-2._wp*GcollL1D(i))
          
          ! !! !Kappa EVDF with spont. and ind. emissions, including EB and CD! !!
          ! ILk1D(i)= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m)*(kappa-1.5_wp)) &
          !         * gamma(kappa)/gamma(kappa-0.5_wp) &
          !         * (1._wp+Zlq2/Q2/VRTeTs(m)/(Kappa-1.5_wp))**(-kappa)+Q*Q2*BremL1D(i)) &
          !         / (2._wp*sqrt(Pi)*VRNeNS(m)/(VRTeTs(m)*(kappa-1.5_wp))**(1.5_wp)*(Zlq2) &
          !         * gamma(kappa+1._wp)/gamma(kappa-0.5_wp) &
          !         * (1._wp+Zlq2/Q2/VRTeTs(m)/(kappa-1.5_wp))**(-(kappa+1._wp)) &
          !         - 2._wp*Q*Q2*GcollL1D(i))
          
          !! !Maxwellian EVDF with spont. and ind. emissions, including EB and CD! !!
          IL1D(i)= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m)) &
               * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q + BremL1D(i)) &
               / (2._wp*sqrt(Pi)*VRNeNS(m)/sqrt(VRTeTs(m)**3)*(Zlq2) &
               * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q - 2._wp*GcollL1D(i))
       end do
       Q= 0._wp
       Q2= 0._wp
       Zlq= 0._wp
       Zlq2= 0._wp
       do i=1,nqx
          Qx=VQx(i)
          do j=1,nqz
             Qz=VQz(j)
             Q2=Qx**2+Qz**2
             Q=sqrt(Q2)
             if(Q<=5.e-3_wp) Q=5.e-3_wp
             call Locate(VQQ,nqcd,Q,ires)
             call Aitp1d2(nqcd,VQQ,Vint1D,Q,Aux,ires)
             Vint(i,j)=Aux
             call Locate(VQQ,nqcd,Q,ires)
             call Aitp1d2(nqcd,VQQ,IL1D,Q,Aux,ires)
             IL(i,j)=Aux
             ! call Locate(VQQ,nqcd,Q,ires)
             ! call Aitp1d2(nqcd,VQQ,ILk1D,Q,Aux,ires)
             ! ILk(i,j)=Aux
             ! call Locate(VQQ,nqcd,Q,ires)
             ! call Aitp1d2(nqcd,VQQ,ILMaxw1D,Q,Aux,ires)
             ! ILMaxw(i,j)=Aux
             Zlq= ZL(Qx,Qz)
             Zlq2= Zlq**2
             !! !Initial spectrum for Maxwellian EVDF w/ spont. and ind. emmissions! !!
             IL0(i,j)=Geff*VRNeNs(m)*VRTeTs(m)/2._wp/Zlq2
             
             !! !Initial spectrum for kappa EVDF, w/ brem and gcoll! !!
             ! IL(i,j)= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m)*(Kappa-1.5_wp)) &
             !      * gamma(Kappa)/gamma(Kappa-0.5_wp) &
             !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(Kappa-1.5_wp))**(-Kappa)+Q*Q2*BremLp(i,j)) &
             !      / (2._wp*sqrt(Pi)*VRNeNS(m)/(VRTeTs(m)*(Kappa-1.5_wp))**(1.5_wp)*(Zlq2) &
             !      * gamma(Kappa+1._wp)/gamma(Kappa-0.5_wp) &
             !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(Kappa-1.5_wp))**(-(Kappa+1._wp)) &
             !      - 2._wp*Q*Q2*GcollLp(i,j))
             
             !! !Initial spectrum for Maxwellian EVDF, w/ brem and gcoll! !!
             ! IL(i,j)= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m)) &
             !      * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q + BremLp(i,j)) &
             !      / (2._wp*sqrt(Pi)*VRNeNS(m)/sqrt(VRTeTs(m)**3)*(Zlq2) &
             !      * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q - 2._wp*GcollLp(i,j))
             ! BrGcL(i,j)= -BremLp(i,j)/GcollLp(i,j)
          end do
       end do
       call output("IL ")
       !call output("ILM")
       call output("GcL")
       call output("BrL")
       call output("Vit")
    end select
    return
  end subroutine wave_init

  real(wp) function ZL(Qx,Qz)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    implicit none
    real(wp), intent(in) :: Qx,Qz
    real(wp) :: Q2
    integer :: m

    m= 1
    Q2= Qx**2+Qz**2
    ZL= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
    if (Qz < 0._wp) then
       ZL= -ZL
    else
    end if
    return
  end function ZL
  
  subroutine coll_damping
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    use math_constants
    use phys_constants
    use DQsimp_m
    use zetaFunc_m
    implicit none
    real(wp) :: Q,Q2
    real(wp) :: Qx,Qz
    real(wp) :: Zlq,Zlq2
    real(wp) :: AuxSig,Aux
    real(wp) :: Res,Res1L,Res2L
    integer :: i,j,m
    integer :: ires,sigpm
    
    m= 1
    Geff= G/(2._wp*sqrt(2._wp)*(4._wp*Pi)**2)

    ! do i=1,nqcd
    !    Q= VQQ(i)
    !    Q2= Q*Q
    !    Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
    !    Zlq2= Zlq**2
    !    Aux1_Gcoll(1)= Q
    !    Aux1_Gcoll(2)= Q2
    !    Aux1_Gcoll(3)= Zlq
    !    Aux1_Gcoll(4)= Zlq2
    !    call dqsimp(Aux_GcollL,1.e-4_wp,4._wp,res1L)
    !    call dqsimpb(Aux_GcollL,4._wp,1.e+30_wp,res2L)
    !    res= Res1L+Res2L
    !    !! !"New" approximation! !!
    !    ! GcollL1D(i)= -8._wp*sqrt(Pi**3)*VRNeNs(m)**4/sqrt(VRTeTs(m)**7) &
    !    !      * Geff/Zlq**2/Q2/Q2*Res
    !    GcollL1D(i)= -8._wp*sqrt(Pi)*VRNeNs(m)**4/sqrt(VRTeTs(m)**7) &
    !         * Geff/Zlq**2/Q2*Res
    ! end do

    do i=1,nqcd
       Q= VQQ(i)
       Q2= Q**2
       Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
       Zlq2= Zlq**2
       GcollL1D(i)= -4._wp*sqrt(Pi)*VRNeNs(m)**4/sqrt(VRTeTs(m)**7) &
            * Geff/Zlq**2/Q2*Vint1D(i)
       GcollLp1D(i)= -16._wp*sqrt(Pi)*VRNeNs(m)**4/sqrt(VRTeTs(m)**7) &
            * Geff*Zlq**2/Q2*Vint1D(i)
    end do
    
    do i= 1,nqx
       Qx= VQx(i)
       do j= 1,nqz
          Qz= VQz(j)
          Q2= Qx**2+Qz**2
          Q= sqrt(Q2)
          if(Q<=5.e-3_wp) Q=5.e-3_wp
          call Locate(VQQ,nqcd,Q,ires)
          call Aitp1d2(nqcd,VQQ,GcollL1D,Q,Aux,ires)
          GcollLp(i,j)= Aux
          GcollLm(i,j)= Aux
       end do
    end do
    return
  end subroutine coll_damping

  ! real(wp) function Aux_GcollL(Qp)
  !   use,intrinsic:: iso_fortran_env, only: wp=>real64
  !   use Common_Params
  !   use Common_Arrays
  !   use Math_Constants
  !   use DQsimp_m
  !   use zetaFunc_m
  !   implicit none
  !   real(wp), intent(in) :: Qp
  !   real(wp) :: Q,Q2,Qp2
  !   real(wp) :: Zlq,Zlq2,Aux
  !   real(wp) :: AuxA,AuxB,AuxA2,AuxB2
  !   real(wp) :: AuxC,AuxC2
  !   real(wp) :: Qp2EpsQpSq
  !   complex(wp) :: Qp2EpsQp
  !   complex(wp) :: Zeta
  !   integer :: m
    
  !   m= 1 
  !   Q= Aux1_Gcoll(1)
  !   Q2= Aux1_Gcoll(2)
  !   Zlq= Aux1_Gcoll(3)
  !   Zlq2= Aux1_Gcoll(4)
  !   ! if(Q > 2.e-4_wp) then
  !   !    if(Qp > 2.e-4_wp) then
  !   Qp2= Qp**2
  !   Zeta= Zlq/Qp/sqrt(VRTeTs(m))
  !   Qp2EpsQp= Qp2+2._wp*VRNeNs(m)/VRTeTs(m) *(1._wp+Zeta*Zfn(Zeta))
  !   Qp2EpsQpSq= (abs(Qp2EpsQp))**2
  !   AuxA= -2._wp*Q*Qp
  !   AuxB= 2._wp*(1._wp+RTeTi)*VRNeNs(m)/VRTeTs(m)+Q2+Qp2
  !   AuxA2= AuxA**2
  !   AuxB2= AuxB**2
  !   Aux= (Qp*(-AuxA2+2._wp*AuxB2)/(AuxB2-AuxA2)+(AuxB/2._wp/Q) &
  !        * log(1._wp+2._wp*AuxA/(AuxB-AuxA)))
  !   Aux_GcollL= exp(-Zlq2/Qp2/VRTeTs(m))*Aux/Qp2EpsQpSq
  !   ! else
  !   !       Aux_GcollL= 0._wp
  !   !    end if
  !   ! else
  !   !    Aux_GcollL= 0._wp
  !   ! end if
  !   return
  ! end function Aux_GcollL

  subroutine Bremsstrahlung
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    use DQsimp_m
    implicit none
    real(wp) :: Qx,Qz
    real(wp) :: Q,Q2
    real(wp) :: Zlq,Zlq2
    real(wp) :: Aux,ResL
    real(wp) :: Res1L,Res2L
    integer :: i,j
    integer :: m,ires
    
    m= 1
    Geff= G/(2._wp*sqrt(2._wp)*(4._wp*Pi)**2)
    
    ! BremL1D= 0._wp
    ! do i= 1,nqcd
    !    Q= VQQ(i)
    !    Q2= Q**2
    !    Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
    !    Zlq2=Zlq**2
    !    Aux1_Bremss(1)= Q
    !    Aux1_Bremss(2)= Q2
    !    Aux1_Bremss(3)= Zlq
    !    Aux1_Bremss(4)= Zlq2

    !    !! !Starting from this loop we have the previous approximation! !!
    !    ! do j=1,nmu
    !    !    Mu= Vmu(j)
    !    !    Aux1_Bremss(5)= Mu
    !    ! call DQsimp(Aux_BremL,1.e-4_wp,4._wp,Res1L)
    !    ! call DQsimpb(Aux_BremL,4._wp,1.e+30_wp,Res2L)
    !    ! VintL(j)=Res1L+Res2L
    !    ! end do
    !    ! call Simpson(Vmu,VintL,nmu,ResL)
    !    ! BremL1D(i)=384._wp*sqrt(pi)/Zlq2*(1._wp-1._wp/RMiMe/VRTiTs(m))**2 &
    !    !      * VRNeNs(m)**4/VRTeTs(m)*Geff**2/Q2*ResL
    !    !! !End of the previous approximation! !!
       
    !    !! !"New" approximation! !!
    !    call DQsimp(Aux_BremL,1.e-4_wp,4._wp,Res1L)
    !    call DQsimpb(Aux_BremL,4._wp,1.e-30_wp,Res2L)
    !    ResL=Res1L+Res2L
    !    BremL1D(i)= 24.0_wp*sqrt(Pi)*VRNeNs(m)**4/sqrt(VRTeTs(m)**5) &
    !         * Geff**2/Zlq**2/Q2*ResL
    ! end do

    do i=1,nqcd
       Q= VQQ(i)
       Q2= Q**2
       Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
       Zlq2= Zlq**2
       BremL1D(i)= 24.0_wp*sqrt(Pi)*VRNeNs(m)**4/sqrt(VRTeTs(m)**5) &
            * Geff**2/Zlq**2/Q2*Vint1D(i)
    end do
    
    do i= 1,nqx
       Qx= VQx(i)
       do j= 1,nqz
          Qz= VQz(j)
          Q2= Qx**2+Qz**2
          Q= sqrt(Q2)
          if(Q<=5.e-3_wp) Q=5.e-3_wp
          call Locate(VQQ,nqcd,Q,ires)
          call Aitp1d2(nqcd,VQQ,BremL1D,Q,Aux,ires)
          BremLp(i,j)= Aux
          BremLm(i,j)= Aux
       end do
    end do
    return
  end subroutine Bremsstrahlung
  
  ! real(wp) function Aux_BremL(Qp)
  !   use,intrinsic:: iso_fortran_env, only: wp=>real64
  !   use Common_Params
  !   use Common_Arrays
  !   use Math_Constants
  !   ! use DQsimp_m
  !   use zetaFunc_m
  !   implicit none
  !   real(wp), intent(in) :: Qp
  !   real(wp) :: Q,Q2
  !   real(wp) :: Qp2,QQpMu
  !   real(wp) :: Qp2EpsQpSq
  !   real(wp) :: Zlq,Zlq2,Aux
  !   real(wp) :: AuxA,AuxB,AuxA2,AuxB2
  !   real(wp) :: Mu,Aux0
  !   complex(wp) :: Qp2EpsQp
  !   complex(wp) :: Zeta
  !   integer :: m

  !   m= 1
  !   Aux0= 0._wp
  !   Aux_BremL=0._wp
    
  !   Q= Aux1_Bremss(1)
  !   Q2= Aux1_Bremss(2)
  !   Zlq= Aux1_Bremss(3)
  !   Zlq2=Aux1_Bremss(4)

  !   !! !Previous approximation! !!
  !   ! Mu= Aux1_Bremss(5)
  !   ! if(Q > 2.e-4_wp) then
  !   !    if(Qp > 2.e-4_wp) then
  !   !       Qp2= Qp**2
  !   !       QQpMu=Q*Qp*Mu
  !   !       Aux0= Qp2*Qp2*(Q2+Qp2-2._wp*QQpMu) &
  !   !            / (2._wp/VRTeTs(m)*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Qp2)**2 &
  !   !            / (2._wp/VRTeTs(m)*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Q2+Qp2-2._wp*QQpMu)**2
  !   !       Aux_BremL= Aux0*(sqrt(1._wp/(VRTeTs(m)*(1._wp*Qp2 &   ! b=e e a=e
  !   !            + 1._wp*(Q2+Qp2-2._wp*QQpMu)))) &    
  !   !            * exp(-Zlq2/(VRTeTs(m)*(1._wp*Qp2 &
  !   !            + 1._wp*(Q2+Qp2-2._wp*QQpMu)))) &
  !   !            + sqrt(1._wp/(VRTeTs(m)*(1._wp*Qp2 &   ! b=e e a=i
  !   !            + VRTiTs(m)/RMiMe*(Q2+Qp2-2._wp*QQpMu)))) &
  !   !            * exp(-Zlq2/(VRTeTs(m)*(1._wp*Qp2 &
  !   !            + VRTiTs(m)/RMiMe*(Q2+Qp2-2._wp*QQpMu)))) &
  !   !            + sqrt(1._wp/(VRTeTs(m)*(VRTiTs(m)/RMiMe*Qp2 &  !b=i e a=e
  !   !            + 1._wp*(Q2+Qp2-2._wp*QQpMu)))) &
  !   !            * exp(-Zlq2/(VRTeTs(m)*(VRTiTs(m)/RMiMe*Qp2 &
  !   !            + 1._wp*(Q2+Qp2-2._wp*QQpMu)))) &
  !   !            + sqrt(1._wp/(VRTeTs(m)*(VRTiTs(m)/RMiMe*Qp2 &  !b=i e a=i
  !   !            + VRTiTs(m)/RMiMe*(Q2+Qp2-2._wp*QQpMu)))) &
  !   !            * exp(-Zlq2/(VRTeTs(m)*(VRTiTs(m)/RMiMe*Qp2 &
  !   !            + VRTiTs(m)/RMiMe*(Q2+Qp2-2._wp*QQpMu)))))
  !   !    else
  !   !       Aux_BremL= 0._wp
  !   !    end if
  !   ! else
  !   !    Aux_BremL= 0._wp
  !   ! end if

  !   !! !New approximation! !!
  !   ! if(Q > 2.e-4_wp) then
  !   !    if(Qp > 2.e-4_wp) then
  !   Qp2= Qp**2
  !   Zeta= Zlq/Qp/sqrt(VRTeTs(m))
  !   Qp2EpsQp= Qp2+2._wp*VRNeNs(m)/VRTeTs(m) *(1._wp+Zeta*Zfn(Zeta))
  !   Qp2EpsQpSq= (abs(Qp2EpsQp))**2
  !   AuxA= -2._wp*Q*Qp
  !   AuxB= 2._wp*(1._wp+RTeTi)*VRNeNs(m)/VRTeTs(m)+Q2+Qp2
  !   AuxA2= AuxA**2
  !   AuxB2= AuxB**2
  !   Aux= (Qp*(-AuxA2+2._wp*AuxB2)/(AuxB2-AuxA2)+(AuxB/2._wp/Q) &
  !        * log(1._wp+2._wp*AuxA/(AuxB-AuxA)))
  !   Aux_BremL= exp(-Zlq2/Qp2/VRTeTs(m))*Aux/Qp2EpsQpSq
  !      ! else
  !   !       Aux_BremL= 0._wp
  !   !    end if
  !   ! else
  !   !    Aux_BremL= 0._wp
  !   ! end if
    
  !   return
  ! end function Aux_BremL

  subroutine Simpson(Vx,F,N,Res)
    ! Version Fortran 95, Aug 2006.
    ! The Jun 2018 version is suitable for
    ! both odd and even values of N.
    ! If N is odd, the subroutine uses the
    ! regular Simpson formula.
    ! If N is even, the subroutine uses the
    ! regular Simpson rule for i=2:N-4 and
    ! the 3/8 Simpson formula for i=N-4,N.
    !! use precision_m
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: N,Nm1,Nm2,Nm3,Nm4,Nm5,i
    real(wp), dimension(N) :: F,Vx
    real(wp) :: Res,Res1,H

    if ( N < 5 ) then
       open(1,FILE='Warning_Simpson.wt')
       if(mod(N,2) /= 0) then
          write(1,*) ' N must be larger or equal 5!'
       else
          write(1,*) ' N must be larger or equal 6!'
       end if
       close(1)
       stop
    else
    end if

    H= ( Vx(N) - Vx(1) )/(N-1)

    Nm1= N-1
    Nm2= N-2
    Nm3= N-3
    Nm4= N-4
    Nm5= N-5
    Res= 0._wp
    Res1= 0._wp
    if(mod(N,2) .ne. 0) then
       ! print*, "odd"
       do i= 2,Nm1,2
          Res= Res + 4._wp*F(i)
       end do
       do i= 3,Nm2,2
          Res= Res + 2._wp*F(i)
       end do
       Res= Res + F(1) + F(N)
       Res= Res*H/3._wp
    else
       ! print*, "even"
       do i= 2,Nm4,2
          Res= Res + 4._wp*F(i)
       end do
       do i= 3,Nm5,2
          Res= Res + 2._wp*F(i)
       end do
       Res= Res + F(1) + F(Nm3)
       Res1= F(Nm3)+ 3._wp*(F(Nm2)+F(Nm1))+F(N)
       Res= Res*H/3._wp + Res1*3._wp*H/8._wp
    end if
    return

  end subroutine Simpson
  
  subroutine Aitp1d2(nx,Vx,Fx,Xp,Fp,i)
    ! Uses linear interpolation in order to obtain the value of a function
    ! F, at the point Xp. The function F is given as a set of points Fx(x).
    ! Uses subroutine Locate (Numerical Recipes, P. 96)
    ! Version Fortran 95, Aug, 2006.
    ! Version of Aitp1d, without the call to Locate, which is called outside.
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: nx,i
    real(wp), dimension(nx) :: Vx,Fx
    real(wp) :: Xp,Fp,Aux
    
    !CALL Locate(Vx,nx,Xp,i)
    
    if ( i<=0 .or. i>=nx ) then
       FP= 0._wp
    else
       Aux= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
       Fp= Fx(i) + ( Fx(i+1)-Fx(i) )*Aux
    end if
    return
  end subroutine Aitp1d2
  
  subroutine Locate(Xx,N,X,J)
    ! Given an array Xx of lenght N, and given a value X, returns a value 
    ! J such that X is between Xx(J) and Xx(J+1).
    ! Xx must be monotonic, either increasing or decreasing.
    ! J=0 or J=N is returned to indicate that X is out of range.
    ! See NUMERICAL RECIPES.
    ! Version Fortran 95, Aug 2006.
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: N,J,JL,JU,JM
    real(wp), dimension (N) :: Xx
    real(wp) :: X
    
    JL= 0
    JU= N+1
10  if ( JU-JL .gt. 1 ) then
       JM= ( JU+JL ) / 2
       if ( (Xx(N).gt.Xx(1)).eqv.(X.gt.Xx(JM)) ) then
          JL= JM
       else
          JU= JM
       end if
       GO TO 10
    end if
    J= JL
    return
  end subroutine Locate

  subroutine output(WriteThis)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use math_constants
    use common_params
    use common_arrays
    implicit none
    real(wp) :: Aux,Dqx,Dqz
    real(wp) :: Qx,Qz,Q
    integer :: i,k
    integer, parameter :: Step=1
    character(LEN=3) :: WriteThis
    Geff= G/(2._wp*sqrt(2._wp)*(4._wp*Pi)**2)
    
    select case(WriteThis)
    case("Fe ")
       open(1,FILE='Fe.wt')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe(nux+1-i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe(nux+1-i,k)
             end if
             write(1,*) -VUx(nux+1-i),VUz(k),Aux
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe(i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe(i,k)
             end if
             write(1,*) VUx(i),VUz(k),Aux
          end do
       end do
       close(1)
       open(1,FILE= "Fe1D.wt")
       do k= 1,nuz
          if(abs(Fe(1,k))<=EpsMin) then
             Aux= EpsMin
          else
             Aux= Fe(1,k)
          end if
          write(1,*) VUz(k),Aux
       end do
       close(1)
       open(1,FILE='Fe0.wt')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe0(nux+1-i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe0(nux+1-i,k)
             end if
             write(1,*) -VUx(nux+1-i),VUz(k),Aux
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe0(i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe0(i,k)
             end if
             write(1,*) VUx(i),VUz(k),Aux
          end do
       end do
       close(1)
       open(1,FILE= "Fe01D.wt")
       do k= 1,nuz
          if(abs(Fe0(1,k))<=EpsMin) then
             Aux= EpsMin
          else
             Aux= Fe0(1,k)
          end if
          write(1,*) VUz(k),Aux
       end do
       close(1)
       ! open(1,FILE='Fek.wt')
       ! do i= 1,nux,Step
       !    write(1,*)' '
       !    do k= 1,nuz,Step
       !       if(abs(Fek(nux+1-i,k))<=EpsMin) then
       !          Aux= EpsMin
       !       else
       !          Aux= Fek(nux+1-i,k)
       !       end if
       !       write(1,*) -VUx(nux+1-i),VUz(k),Aux
       !    end do
       ! end do
       ! do i= 1,nux,Step
       !    write(1,*)' '
       !    do k= 1,nuz,Step
       !       if(abs(Fek(i,k))<=EpsMin) then
       !          Aux= EpsMin
       !       else
       !          Aux= Fek(i,k)
       !       end if
       !       write(1,*) VUx(i),VUz(k),Aux
       !    end do
       ! end do
       ! close(1)
       ! open(1,FILE= "Fek1D.wt")
       ! do k= 1,nuz
       !    if(abs(Fek(1,k))<=EpsMin) then
       !       Aux= EpsMin
       !    else
       !       Aux= Fek(1,k)
       !    end if
       !    write(1,*) VUz(k),Aux
       ! end do
       ! close(1)
    case("IL ")
       open(1,FILE='IL0.wt')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),IL0(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),VQz(k),IL0(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) VQx(i),-VQz(nqz+1-k),IL0(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) VQx(i),VQz(k),IL0(i,k)
          end do
       end do
       close(1)
       open(1,FILE= "IL01D.wt")
       do k= 1,nqcd
          write(1,*) VQQ(k),IL01D(k)
       end do
       close(1)
       open(1,FILE='IL.wt')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),IL(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),VQz(k),IL(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) VQx(i),-VQz(nqz+1-k),IL(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) VQx(i),VQz(k),IL(i,k)
          end do
       end do
       close(1)
       open(1,FILE= "IL1D.wt")
       do k= 1,nqcd
          write(1,*) VQQ(k),IL1D(k)
       end do
       close(1)
       ! open(1,FILE='ILk.wt')
       ! do i= 1,nqx,Step
       !    write(1,*)' '
       !    do k= 1,nqz,Step
       !       write(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),ILk(nqx+1-i,nqz+1-k)
       !    end do
       !    do k= 1,nqz,Step
       !       write(1,*) -VQx(nqx+1-i),VQz(k),ILk(nqx+1-i,k)
       !    end do
       ! end do
       ! do i= 1,nqx,Step
       !    write(1,*)' '
       !    do k= 1,nqz,Step
       !       write(1,*) VQx(i),-VQz(nqz+1-k),ILk(i,nqz+1-k)
       !    end do
       !    do k= 1,nqz,Step
       !       write(1,*) VQx(i),VQz(k),ILk(i,k)
       !    end do
       ! end do
       ! close(1)
       ! open(1,FILE= "ILk1D.wt")
       ! do k= 1,nqcd
       !    write(1,*) VQQ(k),ILk1D(k)
       ! end do
       ! close(1)
    case("ILM")
       open(1,FILE='IL_BrGc-Maxw.wt')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),ILMaxw(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),VQz(k),ILMaxw(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) VQx(i),-VQz(nqz+1-k),ILMaxw(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) VQx(i),VQz(k),ILMaxw(i,k)
          end do
       end do
       close(1)
       open(1,FILE= "IL1D_BrGc-Maxw.wt")
       do k= 1,nqcd
          write(1,*) VQQ(k),ILMaxw1D(k)
       end do
       close(1)
    case("GcL")
       open(1,FILE='GcollL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,500) -Qx,-Qz,GcollLm(nqx+1-i,nqz+1-k), &
                  GcollLm(nqx+1-i,nqz+1-k)/Geff
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,500) -Qx,Qz,GcollLp(nqx+1-i,k), &
                  GcollLp(nqx+1-i,k)/Geff
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,500) Qx,-Qz,GcollLm(i,nqz+1-k), &
                  GcollLm(i,nqz+1-k)/Geff
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,500) Qx,Qz,GcollLp(i,k), &
                  GcollLp(i,k)/Geff
          end do
       end do
       close(1)
       open(1,FILE='GcollL1D.wt')
       do i= 1,nqcd
          write(1,500) VQQ(i),GcollL1D(i),GcollL1D(i)/Geff,GcollL1D(i)/GcollLp1D(i)
       end do
       close(1)
       open(1,FILE='GcollLp1D.wt')
       do i= 1,nqcd
          write(1,500) VQQ(i),GcollLp1D(i),GcollLp1D(i)/Geff
       end do
       close(1)
    case("BrL")
       open(1,FILE='BremL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,500) -Qx,-Qz,BremLm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,500) -Qx,Qz,BremLp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,500) Qx,-Qz,BremLm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,500) Qx,Qz,BremLp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='BremL1D.wt')
       do i= 1,nqcd
          write(1,500) VQQ(i),BremL1D(i)
       end do
       close(1)
       case("Vit")
       open(1,FILE='VintL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,Vint(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,Vint(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,Vint(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,Vint(i,k)
          end do
       end do
       close(1)
       open(1,FILE='VintL1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),Vint1D(i)
       end do
       close(1)
    end select
500 format(1x,5(e13.5e3,1x))
    return
  end subroutine output
end module sub_prog
   
program new_eff
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  use common_params
  use common_arrays
  use math_constants
  use phys_constants
  use sub_prog
  implicit none

  ! open(1,FILE='init.wt')
  ! read(1,*) Kappa
  ! read(1,*) H
  ! read(1,*) RTeTi 
  ! read(1,*) G
  ! read(1,*) Ve2C2
  ! read(1,*) BremGcollL
  call allocate_arrays
  call init_vec
  call definitions
  call fe_init
  call wave_init
  ! close(1)
  
end program new_eff
