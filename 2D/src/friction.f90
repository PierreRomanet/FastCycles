module friction
use variables, only: Dc,a,b,f0,V0,fw,Vw,n_fh,nb_element,dp
contains




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Friction Coefficient: Regularised Rate and state law (Lapusta et al., 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rate_state_R(f,f_dot,V_i,theta_i,nb)
integer										:: nb
real(kind=dp), dimension(nb_element)        :: f,f_dot
real(kind=dp), dimension(nb_element)        :: V_i
real(kind=dp), dimension(nb_element)        :: theta_i
! Internal variable
real(kind=dp), dimension(nb_element)        :: coeff

coeff = exp(f0/a+(b/a)*log(V0*theta_e/Dc))
f = a*asinh(V_e/(2._dp*V0)*coeff1)
f_dot = a/sqrt(1._dp/coeff1**2+(V_e/(2*V0))**2)/(2._dp*V0)
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Friction Coefficient: Rate and state law with flash heating (Dunham et al.,2010a)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rate_state_fh(f,f_dot,V_i,theta_i,nb)
integer										:: nb
real(kind=dp), dimension(nb_element)        :: f,f_dot
real(kind=dp), dimension(nb_element)        :: V_i
real(kind=dp), dimension(nb_element)        :: theta_i
! Internal variable
real(kind=dp), dimension(nb_element)        :: fss


fss = f0 + (a-b)*log(V/V0)
fss = fw + (fss-fw)/(1._dp+(V/Vw)**n_fh)**(1._dp/n_fh)
f = fss + b*log(V*theta/Dc)

f_dot = b/V+(a-b)/V/(1._dp+(V/Vw)**n_fh)**(1._dp/n_fh) & 
      -(f0-fw+(a-b)*log(V/V0))*V**(n_fh-1._dp)/Vw**n_fh/(1._dp+(V/Vw)**n_fh)**(1._dp/n_fh+1._dp)

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Friction Coefficient: Rate and state law 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rate_state(f,f_dot,V_i,theta_i,nb)
integer										:: nb
real(kind=dp), dimension(nb_element)        :: f,f_dot
real(kind=dp), dimension(nb_element)        :: V_i
real(kind=dp), dimension(nb_element)        :: theta_i 


f = f0+a*log(V/V0)+b*log(theta*V0/Dc)
f_dot = a/V
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   State evolution law: slip law
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slip_law(V_i,theta_i,nb)
integer										:: nb
real(kind=dp), dimension(nb_element)        :: slip_law
real(kind=dp), dimension(nb_element)        :: V_i
real(kind=dp), dimension(nb_element)        :: theta_i

slip_law = -V*theta/Dc*log(abs(V)*theta/Dc)

end function


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   State evolution law: aging law
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function aging_law(V_i,theta_i,nb)
integer										:: nb
real(kind=dp), dimension(nb_element)        :: aging_law
real(kind=dp), dimension(nb_element)        :: V_i
real(kind=dp), dimension(nb_element)        :: theta_i

aging_law = 1._dp-abs(V)*theta/Dc

end function
end module







