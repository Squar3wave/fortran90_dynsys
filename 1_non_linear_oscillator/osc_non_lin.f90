program rkprova

  use carepackage
  use rungekutta4

  IMPLICIT NONE
  
  real(dp) :: t_0, dt, t_max, low_b, hi_b, interval
  integer :: N, i, j, l, sigma_size
  real(dp), allocatable, dimension(:) :: x, u, t, sigma, max_p,max_m, x_red
  character(len=30) :: file_name    
  character(len= 3) :: name_helper
  
  open( 1, file= "data_final/data_final_p.dat")  
  open( 2, file= "data_final/data_final_m.dat")
  open(50, file= "maxes/max_p.dat")
  open(51, file= "maxes/max_m.dat")  

  !----------------------------------------------------------------------------------------
  ! Intial variable definition and array allocation
  
  t_0        = 0.0_dp
  t_max      = 800.0_dp
  N          = 1000
  dt         = (t_max-t_0)/N

  low_b      = -9._dp
  hi_b       = 10.5_dp
  interval   =  0.5_dp
  sigma_size = int((hi_b - low_b)/interval) + 1
  
  allocate(    x(0:N-1))
  allocate(    u(0:N-1))
  allocate(    t(0:N-1))
  allocate(x_red(0:199))
  allocate(sigma(0:sigma_size-1))
  allocate(max_p(0:sigma_size-1))
  allocate(max_m(0:sigma_size-1))

  !----------------------------------------------------------------------------------------
  ! Time and sigma values generation using a custom arange rewrite from python  
  
  call arange(t_0, t_max, dt, t)
  call arange(low_b, hi_b, interval, sigma)

  !----------------------------------------------------------------------------------------
  ! Array initialization at 0 for precaution   
  
  x = 0.0_dp
  u = 0.0_dp  

  max_p = 0.0_dp
  max_m = 0.0_dp
  
  !----------------------------------------------------------------------------------------
  ! 1st evolution with Rungew-Kutta 4th order, with crescent values of sigma,
  ! each starting from the last value of the previous,
  ! starting from the second iteration
  
  rk4_crescent_order:do i=0, sigma_size-1
    
    max_p(i) = 0.0_dp

    !--------------------------------------------------------------------------------------
    ! Necessary code to name the files properly  
    
    write(name_helper, '(f0.2)') abs(sigma(i) - int(sigma(i)))
    if(sigma(i)<0) then
      write (file_name,"('data_p/data_p_sigma_',i0.2,A,'.dat')") int(sigma(i)), name_helper  
    else
      write (file_name,"('data_p/data_p_sigma_+',i0.2,A,'.dat')") int(sigma(i)), name_helper   
    endif
    
    open(i+3, file = trim(file_name))

    !--------------------------------------------------------------------------------------
    ! R-K evolution for each value of sigma
    
    call rk4_2ndorder(x(0), u(0), sigma(i), dt, N, t, x, u) 
    
    !--------------------------------------------------------------------------------------
    ! Storing of the last 200 elements of the position, where perturbative
    ! phenomena have subsided, in order to search for maximums
    
    do j= 0, 199
    
      x_red(j) = x(j+800)    
    
    end do
    
    call array_max_finder(x_red,max_p(i))
    write(50,*) sigma(i), max_p(i)
    
    write_on_file_p:do j=0, N-1
    
      write(i+3,*) t(j), x(j)
            
      if (i == sigma_size-1) then
        write(1,*) t(j), x(j)
      end if
    
    end do write_on_file_p

    x(0) = x(N-1)
    u(0) = u(N-1)
    
    close(i+3)
  
  end do rk4_crescent_order

  !----------------------------------------------------------------------------------------
  ! Reinitialization of the needed arrays
  
  x = 0.0_dp
  u = 0.0_dp  
 
  !----------------------------------------------------------------------------------------
  ! 2nd evolution with Rungew-Kutta 4th order, with decrescent values of sigma,
  ! each starting from the last value of the previous,
  ! starting from the second iteration
  
  
  rk4_decrescent_order:do i=0, sigma_size-1
    
    max_m(i) = 0.0_dp
    l = sigma_size-1-i

    !--------------------------------------------------------------------------------------
    ! Necessary code to name the files properly  
    
    write(name_helper, '(f0.2)') abs(sigma(l) - int(sigma(l)))
    if(sigma(l)<0) then
      write (file_name,"('data_m/data_m_sigma_',i0.2,A,'.dat')") int(sigma(l)), name_helper  
    else
      write (file_name,"('data_m/data_m_sigma_+',i0.2,A,'.dat')") int(sigma(l)), name_helper   
    endif
    
    open(i+3, file = trim(file_name))

    !--------------------------------------------------------------------------------------
    ! R-K evolution for each value of sigma    
    
    call rk4_2ndorder(x(0), u(0), sigma(l), dt, N, t, x, u) 

    !--------------------------------------------------------------------------------------
    ! Storing of the last 200 elements of the position, where perturbative
    ! phenomena have subsided, in order to search for maximums
    
    do j= 0, 199
    
      x_red(j) = x(j+800)    
    
    end do
    
    call array_max_finder(x_red,max_m(i))
    write(51,*) sigma(l), max_m(i)
    
    write_on_file_m:do j=0, N-1
    
      write(i+3,*) t(j), x(j)
            
      if (i == sigma_size-1) then
        write(2,*) t(j), x(j)
      end if
    
    end do write_on_file_m

    x(0) = x(N-1)
    u(0) = u(N-1)
    
    close(i+3)
  
  end do rk4_decrescent_order
  
  close( 1)
  close( 2)
  close(50)
  close(51)
  
end program
