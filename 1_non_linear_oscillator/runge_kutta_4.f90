module rungekutta4

use carepackage

implicit none


contains

real(dp) function F(x, v, t, sigma_temp)

  use carepackage
  implicit none      
  real(dp), intent(in) :: x, v, t, sigma_temp 
  real(dp) :: lambda, omega_0, eps, g, ff
  
  eps = 0.1_dp
  omega_0 = 1.0_dp
  g = 2.0_dp
  ff = 10.0_dp  
  lambda = 0.3_dp

  F = -(2*lambda*eps*v + (omega_0**2)*x + eps*g*(x**3)) + eps*ff*cos((omega_0 + eps*sigma_temp)*t)
  
end function F


subroutine rk4_1storder

end subroutine

subroutine rk4_2ndorder (x0_0, x1_0, sigma, h, arrsize, time, x_out, u_out)

  integer,  intent(in)                 :: arrsize
  real(dp), intent(in)                 :: x0_0, x1_0, h, sigma
  real(dp), dimension (:), intent(in)  :: time(0:1000-1)
  real(dp), intent(out)                :: x_out(0:), u_out(0:)
  real(dp), dimension(:)               :: m(0:4-1),k(0:4-1)   
  integer                              :: i
  real(dp), dimension(:) :: vec_coeff(0:3)  
  
  vec_coeff(0) = 0.166666666666666_dp
  vec_coeff(1) = 0.333333333333333_dp
  vec_coeff(2) = 0.333333333333333_dp
  vec_coeff(3) = 0.166666666666666_dp

  x_out(0) = x0_0
  u_out(0) = x1_0
  
  do i=1, arrsize-1
  
    m(0) = h*u_out(i-1)
    
    k(0) = h*F(x_out(i-1), &
               u_out(i-1), &
               time(i), sigma)

               
    m(1) = h*(u_out(i-1)+ 0.5_dp*k(0))
    
    k(1) = h*F(x_out(i-1)+0.5_dp*m(1), &
               u_out(i-1)+0.5_dp*k(0), &
               time(i)+0.5_dp*h, sigma)

               
    m(2) = h*(u_out(i-1) + 0.5_dp*k(1))
    
    k(2) = h*F(x_out(i-1)+0.5_dp*m(2), &
               u_out(i-1)+0.5_dp*k(1), &
               time(i)+0.5_dp*h, sigma)
    
    
    m(3) = h*(u_out(i-1)+ k(2))
    
    k(3) = h*F(x_out(i-1)+m(3), &
               u_out(i-1)+k(2), &
               time(i)+h, sigma)
    

    x_out(i) = x_out(i-1) + dot_product(vec_coeff,m)
    u_out(i) = u_out(i-1) + dot_product(vec_coeff,k)

  end do
 
end subroutine


end module
