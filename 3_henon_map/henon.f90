program henon_map

  use carepackage
  implicit none

  real(dp), dimension(:,:), allocatable :: x, dx, &
                                           x_1, dx_1, &
                                           x_2, dx_2
  real(dp), dimension(:),   allocatable :: pr_12, &
                                           dist, nu, diff, &
                                           ln_dist, ln_nu, c_lin
  integer                               :: N, M, O, i, j ,k
  real(dp)                              :: R, lyap, q, d, d_0, &
                                           R_1, q_1, d_1, lyap_1, &
                                           R_2, q_2, d_2, lyap_2, &
                                           distance, C, FD, FD_C

  !-------------------------------------------------------------------                                           
  ! Variable definitions, array inizialization
  
  N = 10000
  M =  2000
  O =  1000
  
  allocate(    x  (0:1,0:N-1))
  allocate(    x_1(0:1,0:M-1))
  allocate(    x_2(0:1,0:M-1))
  allocate(   dx  (0:1,0:N-1))
  allocate(   dx_1(0:1,0:M-1))
  allocate(   dx_2(0:1,0:M-1))
  allocate(  pr_12(0:1))
  allocate(   dist(0:10))
  allocate(     nu(0:10))
  allocate(ln_dist(0:10))
  allocate(  ln_nu(0:10))
  allocate(   diff(0:O-1))
  allocate(  c_lin(0:O-1))  

  open(1, file ="data/euler_henon_1.dat")
  open(2, file ="data/gp_henon_log.dat")
  open(3, file ="data/gp_henon_coeff.dat")
  
   x     = 0.0_dp
   x_1   = 0.0_dp
   x_2   = 0.0_dp   
  dx     = 0.0_dp
  dx_1   = 0.0_dp
  dx_2   = 0.0_dp
  
  q      = 0.0_dp
  q_1    = 0.0_dp
  q_2    = 0.0_dp
 
  R      = 0.0_dp
  R_1    = 0.0_dp
  R_2    = 0.0_dp

  d      = 0.0_dp
  d_1    = 0.0_dp
  d_2    = 0.0_dp

  d_0    = 0.0_dp 
  
  lyap   = 0.0_dp
  lyap_1 = 0.0_dp
  lyap_2 = 0.0_dp
  
  pr_12  = 0.0_dp
  
  !-------------------------------------------------------------------
  ! Simple evolution of the system via Euler's method
  
   x(0,0) = 0.1_dp
   x(1,0) = 0.1_dp
  
  write(1,*) x(0,0), x(1,0)
  
  do i = 1, N-1
  
     x(:,i) =  henon(x(:,i-1))

    write(1,*) x(0,i), x(1,i)
    
  end do
  
  !-------------------------------------------------------------------
  ! Calculation of Lyapunov exponent with Benettin's alogrithm, 
  ! once again relying on Euler's method for system evolution

  print "('------------------------------------------------------- ')"
  print "('Lyapunov coefficient for single trajectory')"
  print "('calculated with Benettin algorithm')"
  print "('------------------------------------------------------- ')"
  print "(' ')"
  
   x(0,0)  = 0.1_dp
   x(1,0)  = 0.1_dp
  dx(0,0)  = 0.9_dp
  dx(1,0)  = 0.9_dp
  
  d_0      = dot_product(x(:,0),x(:,0))
  
  do i = 1, N-1

     x(:,i) =  henon(x(:,i-1))
    dx(:,i) = dhenon(x(:,i-1),dx(:,i-1))

     if (modulo(i,100) == 0) then
     
       d = dot_product(dx(:,i),dx(:,i))
       q = d/d_0
       R = R + log(sqrt(q))
 
      dx = dx/sqrt(q)
       
     end if

  end do  

  lyap = R/real(N)

  print "('Theoretical value            = ', f11.7)", 0.4191919
  print "('Experimental value           = ', f11.7)", lyap
  
  !-------------------------------------------------------------------
  ! Calculation of Lyapunov exponent for two vectors with Benettin, 
  ! once again relying on Euler's method for system evolution
  ! and orthogonalizing with Gram-Schmidt every n steps 
  ! as required by Benettin's algorithm
  
  print "(' ')"
  print "('------------------------------------------------------- ')"
  print "('Lyapunov coefficients for 2 trajectories ')"
  print "('calculated with Benettin algorithm')"
  print "('------------------------------------------------------- ')"
  print "(' ')"

  
   x_1(0,0) = 0.8_dp
   x_1(1,0) = 0.1_dp
  dx_1(0,0) = 1.0_dp
  dx_1(1,0) = 0.0_dp

   x_2(0,0) = 0.8_dp
   x_2(1,0) = 0.1_dp
  dx_2(0,0) = 0.0_dp
  dx_2(1,0) = 1.0_dp  

  d_0       = dot_product(x_1(:,0),x_1(:,0))
  
  do i = 1, M-1
        
    !-----------------------------------------------------------------
    ! 1st trajectory's i step of the evolution

     x_1(:,i) =  henon(x_1(:,i-1))
    dx_1(:,i) = dhenon(x_1(:,i-1),dx_1(:,i-1))
    
    !-----------------------------------------------------------------
    ! 2nd trajectory's i step of the evolution
    
     x_2(:,i) =  henon(x_2(:,i-1))
    dx_2(:,i) = dhenon(x_2(:,i-1),dx_2(:,i-1))
     
    !-----------------------------------------------------------------
    ! Orthogonalization every 10 steps with subsequent
    ! rinormalization and progressive calculation of
    ! numerators in Lyapunov's formula
    
    if (modulo(i,10) == 0) then

      d_1       = sqrt(dot_product(dx_1(:,i),dx_1(:,i)))
      q_1       = d_1/d_0
      R_1       = R_1 + log(q_1)

      pr_12     = (dot_product(dx_2(:,i),dx_1(:,i))/&
                 dot_product(dx_1(:,i),dx_1(:,i)))*dx_1(:,i)
      
      dx_2(:,i) = dx_2(:,i) - pr_12

      d_2       = sqrt(dot_product(dx_2(:,i),dx_2(:,i)))
      q_2       = d_2/d_0
      R_2       = R_2 + log(q_2)

      dx_1      = dx_1/q_1
      dx_2      = dx_2/q_2
       
    end if

  end do    

  lyap_1 = R_1/real(M)
  lyap_2 = R_2/real(M)

  print "('Lyapunov coefficient for x_1 = ', f11.7)", lyap_1
  print "('Lyapunov coefficient for x_2 = ', f11.7)", lyap_2

  !-------------------------------------------------------------------
  ! Calculation of the fractal dimension using
  ! Grassberg-Procaccia algorithm
  
  print "(' ')"
  print "('------------------------------------------------------- ')"
  print "('Fractal Dimension calculated with Grassberger-Procaccia ')"
  print "('------------------------------------------------------- ')"
  print "(' ')"
  
  dist( 0) = 0.10_dp
  dist( 1) = 0.15_dp
  dist( 2) = 0.20_dp
  dist( 3) = 0.25_dp
  dist( 4) = 0.30_dp
  dist( 5) = 0.35_dp
  dist( 6) = 0.40_dp
  dist( 7) = 0.50_dp
  dist( 8) = 0.60_dp
  dist( 9) = 0.70_dp
  dist(10) = 0.80_dp
  
  do i = 0, 10
  
    C = 0
    
    do j = 0, O-1
      do k = j+1, O-1

        diff     = x(:,j) - x(:,k)
        distance = sqrt(dot_product(diff,diff))  
        
        if (distance < dist(i)) then
        
          C = C+1
        
        end if
        
      end do
    end do
  
    nu(i)      = C/(O**2)
    ln_dist(i) = log(dist(i))
    ln_nu(i)   = log(nu(i))
    
    write(2, *) ln_dist(i), ln_nu(i)
    
  end do
  
  call lin_fit(ln_dist, ln_nu, 11, c_lin)
  
  write(3,*) c_lin(0), c_lin(1)
  
  FD_C = 1.0_dp - (lyap_1/lyap_2)  
  
  FD   = c_lin(1)

  print "('Kaplan-Yorke prediction      = ', f11.7)", FD_C
  print "('Calculated value             = ', f11.7)", FD
  
  close(1)
  close(2)
  close(3)
  
  deallocate(x)
  deallocate(x_1)
  deallocate(x_2)
  deallocate(dx)
  deallocate(dx_1)
  deallocate(dx_2)
  deallocate(pr_12)
  deallocate(dist)
  deallocate(nu)
  deallocate(ln_dist)
  deallocate(ln_nu)
  deallocate(diff)
  deallocate(c_lin)  
  
contains

function henon(x_in)

  real(dp), dimension(:) :: henon(0:1)
  real(dp), intent(in)   :: x_in(0:1)
  real(dp)               :: a,b

  a = 1.4_dp
  b = 0.3_dp

  henon(0) = x_in(1) + 1.0_dp - a*(x_in(0)**2)
  henon(1) = b*x_in(0)

end function henon

function dhenon(x_in, dx_in)

  real(dp)             :: dhenon(0:1)
  real(dp), intent(in) :: x_in(0:1), dx_in(0:1)
  real(dp)             :: a,b

  a = 1.4_dp
  b = 0.3_dp

  dhenon(0) = dx_in(1) - 2*a*x_in(0)*dx_in(0)
  dhenon(1) = b*dx_in(0)

end function

end program
