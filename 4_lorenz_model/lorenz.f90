program lorenz_model

  use carepackage

  IMPLICIT NONE
  
  integer                               :: t_size, i, j, l, ly_c, ly_c_3  
  
  real(dp)                              :: s, r, b,     &
                                           h, t_0, t_f, &
                                           lyap_coeff, &
                                           dist, c, FD, FD_KY
                                           
  real(dp), allocatable, dimension(:)   :: lyap_3, &
                                           dist_l, nu, ln_dis, ln_nu , &
                                           coeff
  
  real(dp), allocatable, dimension(:,:) :: pos_1, dpos_1, &
                                           pos_2, dpos_2, &
                                           pos_a, dpos_a, &
                                           pos_b, dpos_b, &
                                           pos_c, dpos_c, &
                                           stat

  open(1, file = "rk4_lorenz.dat")
  open(2, file = "stz.dat")
  open(3, file = "gp.dat")

  !-----------------------------------------------------------------------------------
  ! legend:
  ! array(i,j) means the first index labels the spatial coordinates,
  ! while the second one labels the time step, if not stated otherwise 
  ! 0==x 1==y 2==z

  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  ! Parameters definition and array allocation

       h = 0.01_dp
     t_0 = 0.0_dp
     t_f = 500.0_dp
  t_size = int((t_f-t_0)/h) + 1

  ly_c   = 10000
  ly_c_3 = 1000000
    
       s = 10.0_dp
       r = 28.0_dp
       b =  8.0_dp/3.0_dp
       
  allocate(lyap_3(0:2))
  allocate( pos_1(0:2,0:t_size-1))
  allocate(dpos_1(0:2,0:t_size-1))
  allocate(  stat(0:2,0:2))
  allocate( pos_2(0:2,0:ly_c-1))
  allocate(dpos_2(0:2,0:ly_c-1))
  allocate( pos_a(0:2,0:ly_c_3-1))
  allocate(dpos_a(0:2,0:ly_c_3-1))
  allocate( pos_b(0:2,0:ly_c_3-1))
  allocate(dpos_b(0:2,0:ly_c_3-1))
  allocate( pos_c(0:2,0:ly_c_3-1))
  allocate(dpos_c(0:2,0:ly_c_3-1))
  allocate(dist_l(0:12))
  allocate(    nu(0:12))  
  allocate(ln_dis(0:12))
  allocate( ln_nu(0:12))
  allocate( coeff(0:1))
  
  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  ! Stationary points calculation
  !
  ! (i,j) i space coord. j stat. point identifier
  
  stat(0,0) = 0.0_dp
  stat(1,0) = 0.0_dp 
  stat(2,0) = 0.0_dp
    
  stat(0,1) = sqrt(b*(r-1.0_dp))
  stat(1,1) = sqrt(b*(r-1.0_dp))  
  stat(2,1) = r-1.0_dp
  
  stat(0,2) = -sqrt(b*(r-1.0_dp))
  stat(1,2) = -sqrt(b*(r-1.0_dp))  
  stat(2,2) = r-1.0_dp

  do i=0,2

    write(2,*) stat(0,i), stat(1,i), stat(2,i)
  
  end do

  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  ! System evolution with 4th order Runge-Kutta
  !
  ! (i,j) i space coord. j time step
  
  pos_1      = 0.0_dp
  
  pos_1(0,0) = 1.0_dp
  pos_1(1,0) = 0.0_dp
  pos_1(2,0) = 0.0_dp

  call rk4_1storder_vec(h, t_size, pos_1)
  
  do i=0,t_size-1
    
    write(1,*) pos_1(0,i), pos_1(1,i), pos_1(2,i)
 
  end do
  
  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  ! Lyapunov exponent calculation with Benettin algorithm,
  ! with new starting conditions
  !
  ! (i,j) i space coord. j time step

  print "('------------------------------------------------------- ')"
  print "('Lyapunov coefficient for single trajectory')"
  print "('calculated with Benettin algorithm')"
  print "('------------------------------------------------------- ')"
  print "(' ')"
  
   pos_2      = 0.0_dp  
  
   pos_2(0,0) = 3.1_dp
   pos_2(1,0) = 4.2_dp
   pos_2(2,0) = 5.7_dp

  dpos_2(0,0) = 1.0_dp
  dpos_2(1,0) = 1.0_dp
  dpos_2(2,0) = 1.0_dp

   lyap_coeff = 0.0_dp
  
  call lyap_calc_rk4(h, ly_c, pos_2, dpos_2, lyap_coeff)

  print "('Lyapunov coefficient for pos_2 = ', f11.7)", lyap_coeff

  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  ! Lyapunov exponent calculation with Benettin algorithm for three trajectories,
  ! orthogonalizing with Gram-Schmidt every n steps as required by the procedure
  
  print "(' ')"
  print "('------------------------------------------------------- ')"
  print "('Lyapunov coefficients for 3 trajectories ')"
  print "('calculated with Benettin algorithm')"
  print "('------------------------------------------------------- ')"
  print "(' ')"

  
   pos_a(0,0) = 1.0_dp
   pos_a(1,0) = 2.0_dp
   pos_a(2,0) = 1.5_dp

  dpos_a(0,0) = 1.0_dp
  dpos_a(1,0) = 0.0_dp
  dpos_a(2,0) = 0.0_dp
  
   pos_b(0,0) = 1.0_dp
   pos_b(1,0) = 2.0_dp
   pos_b(2,0) = 1.5_dp

  dpos_b(0,0) = 0.0_dp
  dpos_b(1,0) = 1.0_dp
  dpos_b(2,0) = 0.0_dp
  
   pos_c(0,0) = 1.0_dp
   pos_c(1,0) = 2.0_dp
   pos_c(2,0) = 1.5_dp

  dpos_c(0,0) = 0.0_dp
  dpos_c(1,0) = 0.0_dp
  dpos_c(2,0) = 1.0_dp
  
            h = 0.001_dp
  
  call lyap_calc_rk4_3 (h, ly_c_3, pos_a, dpos_a, &
                                          pos_b, dpos_b, &
                                          pos_c, dpos_c, &
                                          lyap_3)
                                          
  print "('Lyapunov coefficient for pos_a = ', f11.7)", lyap_3(0)
  print "('Lyapunov coefficient for pos_b = ', f11.7)", lyap_3(1)
  print "('Lyapunov coefficient for pos_c = ', f11.7)", lyap_3(2)
  print "('Calculated sum of the 3 values = ', f11.7)", lyap_3(0) +&
                                                        lyap_3(1) + lyap_3(2)
  print "('Expected sum of the 3 values   = ', f11.7)", -(s+b+1)                                                        

  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  ! Fractal dimension calculation with Grassberger-Procaccia

  print "(' ')"
  print "('------------------------------------------------------- ')"
  print "('Fractal Dimension calculated with Grassberger-Procaccia ')"
  print "('------------------------------------------------------- ')"
  print "(' ')"
  
  dist_l( 0) =  0.2_dp
  dist_l( 1) =  0.4_dp
  dist_l( 2) =  0.7_dp
  dist_l( 3) =  1.0_dp
  dist_l( 4) =  1.5_dp
  dist_l( 5) =  2.0_dp
  dist_l( 6) =  3.0_dp
  dist_l( 7) =  4.0_dp
  dist_l( 8) =  5.0_dp
  dist_l( 9) =  6.0_dp
  dist_l(10) =  7.0_dp
  dist_l(11) =  8.0_dp
  dist_l(12) = 10.0_dp
  
  do i = 0,12
    
    c = 0.0_dp
    
    !$OMP PARALLEL
    
    do j= 0, 1000-1
      do l = j+1, 1000-1
         
        dist = sqrt((pos_a(0,j) - pos_a(0,l))**2 + &
                    (pos_a(1,j) - pos_a(1,l))**2 + &
                    (pos_a(2,j) - pos_a(2,l))**2)
         
        if( dist < dist_l(i) ) then
        
          c = c+1.0_dp
        
        end if
      
      end do  
    end do
    
    !$OMP END PARALLEL

    !---------------------------------------------------------------------------------
    ! Being in presence of a power law we calculate the logarithms for bot x and y
    ! so that we can work with a linear fit
    
    nu(i)     = c/(1000.0_dp**2)
    ln_dis(i) = log(dist_l(i))
    ln_nu(i)  = log(nu(i))

    write(3,*) i, ln_nu(i)
    
  end do
  
  !-----------------------------------------------------------------------------------
  ! Custom linear fit subroutine which returns an array with the coefficients
  ! coeff(0) is the offset, while coeff(1) the angular coefficient
  
  call lin_fit(ln_dis, ln_nu, 13, coeff)
  
  FD_KY = 2.0_dp + (lyap_3(0) + lyap_3(1))/abs(lyap_3(2))
  FD    = coeff(1)

  print "('Theoretical value              = ', f11.7)", 2.08_dp
  print "('Kaplan-Yorke prediction        = ', f11.7)", FD_KY
  print "('Calculated value               = ', f11.7)", FD
  
  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  ! Array deallocation and closing procedures

  close(1)
  close(2)
  close(3)
  
  deallocate(stat)
  deallocate(pos_1)
  deallocate(pos_2)
  deallocate(pos_a)
  deallocate(pos_b)
  deallocate(pos_c)
  deallocate(dpos_1)
  deallocate(dpos_2)
  deallocate(dpos_a)
  deallocate(dpos_b)
  deallocate(dpos_c)
  deallocate(dist_l)
  deallocate(nu)  
  deallocate(ln_dis)
  deallocate(ln_nu)
  deallocate(coeff)

contains

function lorenz (s_in, r_in, b_in,v_in)

  !-----------------------------------------------------------------------------------
  ! 1st order differential equations of the system, 
  ! to be integrated with Runge-Kutta
  
  real(dp), intent(in)                 :: r_in, s_in, b_in
  real(dp), dimension (:), intent(in)  :: v_in(0:2)
  real(dp), dimension(:)               :: lorenz(0:2)
  
  lorenz(0) =  s_in*(-v_in(0) + v_in(1))
  
  lorenz(1) =  r_in*v_in(0) - v_in(1) - v_in(0)*v_in(2)
  
  lorenz(2) = -b_in*v_in(2) + v_in(0)*v_in(1)
 
end function lorenz

subroutine rk4_1storder_vec (h_in, rk_size_in, pos_out)

  !-----------------------------------------------------------------------------------
  ! 4th order Runge-Kutta algorithm calculation

  real(dp), intent(out) :: pos_out(0:,0:)
  real(dp), intent(out) :: h_in
  integer :: i
  integer, intent(in) :: rk_size_in
  real(dp), dimension(:) :: vec_coeff(0:3)  
  
  vec_coeff(0) = 0.166666666666666_dp
  vec_coeff(1) = 0.333333333333333_dp
  vec_coeff(2) = 0.333333333333333_dp
  vec_coeff(3) = 0.166666666666666_dp

  do i=1, rk_size_in-1

    call rk4_calc(vec_coeff, h_in, pos_out(:,i-1), pos_out(:,i))
    
  end do 
 
end subroutine

subroutine rk4_calc(vec_coeff_in, h_in, pos_in, rk4_c)

  real(dp), intent(out)                 :: rk4_c(0:2),vec_coeff_in(0:3)
  real(dp), dimension(:), intent(in)    :: pos_in(0:2)
  real(dp), intent(in)                  :: h_in
  real(dp)                              :: pos_temp(0:2),k(0:2,0:3) 

    !---------------------------------------------------------------------------------
    ! Step 1

    k(:,0)   = h_in * lorenz(s,r,b,pos_in)
    pos_temp = pos_in + k(:,0)*0.5_dp

    !---------------------------------------------------------------------------------
    ! Step 2

    k(:,1)   = h_in * lorenz(s,r,b,pos_temp)
    pos_temp = pos_in + k(:,1)*0.5_dp
    
    !---------------------------------------------------------------------------------
    ! Step 3
   
    k(:,2)   = h_in * lorenz(s,r,b,pos_temp)
    pos_temp = pos_in + k(:,2)
 
    !---------------------------------------------------------------------------------
    ! Step 4
    
    k(:,3)   = h_in * lorenz(s,r,b,pos_temp)
    
    rk4_c(0) = pos_in(0) + dot_product(vec_coeff_in, slicer_k(k,0))
    rk4_c(1) = pos_in(1) + dot_product(vec_coeff_in, slicer_k(k,1))
    rk4_c(2) = pos_in(2) + dot_product(vec_coeff_in, slicer_k(k,2))


end subroutine

subroutine drk4_calc(vec_coeff_in, h_in, pos_in, dpos_in, drk4_c)

  real(dp), intent(out)                 :: drk4_c(0:2),vec_coeff_in(0:3)
  real(dp), dimension(:), intent(in)    :: pos_in(0:2), dpos_in(0:2)
  real(dp), intent(in)                  :: h_in
  real(dp)                              :: dpos_temp(0:2), dk(0:2,0:3)


    !---------------------------------------------------------------------------------
    ! Step 1
    
    dk(:,0)   = h_in * dlorenz(s,r,b,pos_in, dpos_in)
    dpos_temp = dpos_in(:) + dk(:,0)*0.5_dp

    !---------------------------------------------------------------------------------
    ! Step 2

    dk(:,1)   = h_in * dlorenz(s,r,b,pos_in, dpos_temp)    
    dpos_temp = dpos_in(:) + dk(:,1)*0.5_dp      

    !---------------------------------------------------------------------------------
    ! Step 3
    
    dk(:,2)   = h_in * dlorenz(s,r,b,pos_in, dpos_temp)    
    dpos_temp = dpos_in(:) + dk(:,2)    
    
    !---------------------------------------------------------------------------------
    ! Step 4
    
    dk(:,3)   = h_in * dlorenz(s,r,b,pos_in, dpos_temp)    

    drk4_c(0) = dpos_in(0) + dot_product(vec_coeff_in, slicer_k(dk,0))
    drk4_c(1) = dpos_in(1) + dot_product(vec_coeff_in, slicer_k(dk,1))
    drk4_c(2) = dpos_in(2) + dot_product(vec_coeff_in, slicer_k(dk,2))

end subroutine

function dlorenz (s_in, r_in, b_in,v_in, dv_in)

  real(dp), intent(in)                 :: r_in, s_in, b_in
  real(dp), dimension (:), intent(in)  :: dv_in(0:2), v_in(0:2)
  real(dp), dimension(:)               :: dlorenz(0:2)

  dlorenz(0) =  s_in*(-dv_in(0) + dv_in(1))
  
  dlorenz(1) =  r_in*dv_in(0) - dv_in(1) &
                              - dv_in(0)*v_in(2) &
                              - dv_in(2)*v_in(0)
  
  dlorenz(2) = -b_in*dv_in(2) + dv_in(0)*v_in(1) &
                              + dv_in(1)*v_in(0)

end function dlorenz

subroutine lyap_calc_rk4 (h_in, rk_size_in, &
                          pos_out, dpos_out, lyap)

  real(dp), intent(out)  :: pos_out(0:,0:), dpos_out(0:,0:)
  real(dp), intent(out)  :: h_in, lyap
  integer                :: i
  integer, intent(in)    :: rk_size_in
  real(dp)               :: d, d_0, q, r_l
  real(dp), dimension(:) :: vec_coeff(0:3)  
  
  vec_coeff(0) = 0.166666666666666_dp
  vec_coeff(1) = 0.333333333333333_dp
  vec_coeff(2) = 0.333333333333333_dp
  vec_coeff(3) = 0.166666666666666_dp
  
             d = 0.0_dp
             q = 0.0_dp
           r_l = 0.0_dp
           d_0 = sqrt(dot_product(dpos_out(:,0),dpos_out(:,0)))
           
  !-----------------------------------------------------------------------------------
  ! To calculate Lyapunov coefficient with Benettin's algorithm
  ! we have to run R-K algorithms 
  ! for both the lorentz equations and their derivates
  
  call rk4_1storder_vec (h_in, rk_size_in, pos_out) 

  !-----------------------------------------------------------------------------------
  ! Evolution via Runge-Kutta 4th order for the derivates, 
  ! rinormalization every n steps and Lyapunov numerator's calculation
  
  do i=1, rk_size_in-1

    !---------------------------------------------------------------------------------
    ! R-K i step calculation for dpos
  
    call drk4_calc(vec_coeff, h_in, pos_out(:,i-1),dpos_out(:,i-1), dpos_out(:,i))

    !---------------------------------------------------------------------------------
    ! Rinormalization every 100 steps of the evolution,
    ! and progressive calculation of numerator in Lyapunov's formula
                                    
    if (modulo(i,100) == 0) then
    
              d = sqrt(dot_product(dpos_out(:,i),dpos_out(:,i)))
              q = d/d_0
            r_l = r_l + (log(q)/h_in)
       
      dpos_out  = dpos_out/q      
       
    end if
    
  end do
 
  lyap = r_l/real(rk_size_in)
 
end subroutine

subroutine lyap_calc_rk4_3 (h_in, rk_size_in, &
                            pos_a_out, dpos_a_out, &
                            pos_b_out, dpos_b_out, &
                            pos_c_out, dpos_c_out, &
                            lyap_3)

  real(dp), intent(out)  :: pos_a_out(0:,0:), dpos_a_out(0:,0:), &
                            pos_b_out(0:,0:), dpos_b_out(0:,0:), &
                            pos_c_out(0:,0:), dpos_c_out(0:,0:) 
  real(dp), intent(out)  :: h_in
  integer                :: i
  integer, intent(in)    :: rk_size_in
  real(dp), dimension(:) :: pr_ab(0:2), pr_ac(0:2), pr_bc(0:2), &
                            r_l(0:2),lyap_3(0:2), q(0:2)
  real(dp)               :: d_a, d_b, d_c, d_0
  real(dp), dimension(:) :: vec_coeff(0:3)  
  
  vec_coeff(0) = 0.166666666666666_dp
  vec_coeff(1) = 0.333333333333333_dp
  vec_coeff(2) = 0.333333333333333_dp
  vec_coeff(3) = 0.166666666666666_dp
  
           d_a = 0.0_dp
           d_b = 0.0_dp
           d_c = 0.0_dp
             q = 0.0_dp
           r_l = 0.0_dp

           d_0 = sqrt(dot_product(dpos_a_out(:,0),dpos_a_out(:,0)))
           
  !-----------------------------------------------------------------------------------
  ! To calculate Lyapunov coefficient we have to run R-K algorithms 
  ! for both the lorenz equations and their derivatives for all 3 directions
  
  call rk4_1storder_vec (h_in, rk_size_in, pos_a_out) 
  call rk4_1storder_vec (h_in, rk_size_in, pos_b_out)
  call rk4_1storder_vec (h_in, rk_size_in, pos_c_out)  
  
  do i=1, rk_size_in-1

    !---------------------------------------------------------------------------------
    ! R-K i step calculation for dpos_a
     
    call drk4_calc(vec_coeff, h_in, &
                   pos_a_out(:,i-1),dpos_a_out(:,i-1), dpos_a_out(:,i))

    !---------------------------------------------------------------------------------
    ! R-K i step calculation for dpos_b
     
    call drk4_calc(vec_coeff, h_in, &
                   pos_b_out(:,i-1),dpos_b_out(:,i-1), dpos_b_out(:,i))

    !---------------------------------------------------------------------------------
    ! R-K i step calculation for dpos_c
     
    call drk4_calc(vec_coeff, h_in, &
                   pos_c_out(:,i-1),dpos_c_out(:,i-1), dpos_c_out(:,i))
                                      
    !---------------------------------------------------------------------------------
    ! Orthogonalization via Gram-Schmidt,
    ! progressive calculation of numerators in Lyapunov's formula
    ! and vector rinormalization every 100 steps of the evolution    
    
    if (modulo(i,100) == 0) then
                       
                  d_a = sqrt(dot_product(dpos_a_out(:,i),dpos_a_out(:,i)))
                 q(0) = d_a/d_0
               r_l(0) = r_l(0) + (log(q(0))/(h_in))

                pr_ab = (dot_product(dpos_b_out(:,i),dpos_a_out(:,i))/&
                         dot_product(dpos_a_out(:,i),dpos_a_out(:,i)))*dpos_a_out(:,i)

      dpos_b_out(:,i) = dpos_b_out(:,i) - pr_ab
                
                  d_b = sqrt(dot_product(dpos_b_out(:,i),dpos_b_out(:,i)))
                 q(1) = d_b/d_0
               r_l(1) = r_l(1) + (log(q(1))/(h_in))

                pr_ac = (dot_product(dpos_c_out(:,i),dpos_a_out(:,i))/&
                         dot_product(dpos_a_out(:,i),dpos_a_out(:,i)))*dpos_a_out(:,i)
                pr_bc = (dot_product(dpos_c_out(:,i),dpos_b_out(:,i))/&
                         dot_product(dpos_b_out(:,i),dpos_b_out(:,i)))*dpos_b_out(:,i)
                       
      dpos_c_out(:,i) = dpos_c_out(:,i) - pr_ac - pr_bc
              
                  d_c = sqrt(dot_product(dpos_c_out(:,i),dpos_c_out(:,i)))
                 q(2) = d_c/d_0
               r_l(2) = r_l(2) + (log(q(2))/(h_in))

      dpos_a_out(:,i) = dpos_a_out(:,i)/q(0)
      dpos_b_out(:,i) = dpos_b_out(:,i)/q(1)      
      dpos_c_out(:,i) = dpos_c_out(:,i)/q(2)      

    end if
    
  end do

   lyap_3 = r_l/real(rk_size_in)
 
end subroutine

end program
