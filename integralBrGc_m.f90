module commom_var
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp), allocatable :: VRNeNs(:),VRTeTs(:)
  real(wp), dimension(7) :: Aux_int
end module commom_var
module integralBrGc_m
contains
  subroutine integral(nq,VQ,RTeTi,sol)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use commom_var
    use DQsimp_m
    use zetaFunc_m
    implicit none
    real(wp) :: RTeTi
    real(wp) :: VQ(nq)
    real(wp) :: sol(nq)
    real(wp) :: Q,Q2
    real(wp) :: Zlq,Zlq2
    real(wp) :: res1,res2
    integer :: i,m,nrz,nq
    nrz=1
    allocate(VRNeNs(nrz),VRTeTs(nrz))
    m=1
    VRNeNs(m)= 1._wp
    VRTeTs(m)= 1._wp
    Aux_int(5)= VRNeNs(m)
    Aux_int(6)= VRTeTs(m)
    Aux_int(7)= RTeTi
    do i=1,nq
       Q= VQ(i)
       Q2= Q*Q
       Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
       Zlq2= Zlq**2
       Aux_int(1)= Q
       Aux_int(2)= Q2
       Aux_int(3)= Zlq
       Aux_int(4)= Zlq2
       call dqsimp(func_int,1.e-4_wp,4._wp,res1)
       call dqsimpb(func_int,4._wp,1.e+30_wp,res2)
       sol(i)= res1+res2
    end do
    return
  end subroutine integral

  real(wp) function func_int(Qp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use commom_var
    use DQsimp_m
    use zetaFunc_m
    implicit none
    real(wp), intent(in) :: Qp
    real(wp) :: Q,Q2,Qp2
    real(wp) :: Zlq,Zlq2,Aux
    real(wp) :: RTeTi
    real(wp) :: AuxA,AuxB,AuxA2,AuxB2
    real(wp) :: AuxC,AuxC2
    real(wp) :: Qp2EpsQpSq
    complex(wp) :: Qp2EpsQp
    complex(wp) :: Zeta
    integer :: m,nrz
    !nrz=1
    !allocate(VRNeNs(nrz),VRTeTs(nrz))
    m= 1 
    Q= Aux_int(1)
    Q2= Aux_int(2)
    Zlq= Aux_int(3)
    Zlq2= Aux_int(4)
    VRNeNs(m)= Aux_int(5)
    VRTeTs(m)= Aux_int(6)
    RTeTi= Aux_int(7)
    
    Qp2= Qp**2
    Zeta= Zlq/Qp/sqrt(VRTeTs(m))
    Qp2EpsQp= Qp2+2._wp*VRNeNs(m)/VRTeTs(m) *(1._wp+Zeta*Zfn(Zeta))
    Qp2EpsQpSq= (abs(Qp2EpsQp))**2
    AuxA= -2._wp*Q*Qp
    AuxB= 2._wp*(1._wp+RTeTi)*VRNeNs(m)/VRTeTs(m)+Q2+Qp2
    AuxA2= AuxA**2
    AuxB2= AuxB**2
    Aux= (Qp*(-AuxA2+2._wp*AuxB2)/(AuxB2-AuxA2)+(AuxB/2._wp/Q) &
         * log(1._wp+2._wp*AuxA/(AuxB-AuxA)))
    func_int= exp(-Zlq2/Qp2/VRTeTs(m))*Aux/Qp2EpsQpSq
    return
  end function Func_Int

end module integralBrGc_m
