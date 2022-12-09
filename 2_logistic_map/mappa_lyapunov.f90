program mappa_lyap

  use carepackage
  implicit none

  real(dp), allocatable, dimension(:) :: muu, itera, der, ln, numeri, lyapunov, lyap_mu
  real(dp) :: mu_low, mu_high, mu_delta, itera_low, itera_high, itera_delta
  real(dp) :: x_n, x_np1, log_sum
  integer :: mu_size, itera_size, itera_size_2, i, j, k, l
  character(len=37) :: file_name

  !-----------------------------------------------------------------------------------
  ! Variable definition and array allocation
  
  mu_high      = 4.0_dp
  mu_low       = 0.1_dp
  mu_delta     = 0.15_dp
  mu_size      = int((mu_high - mu_low)/mu_delta) + 1
  
  itera_high   = 1000.0_dp
  itera_low    = 2.0_dp
  itera_delta  = 25.0_dp
  itera_size   = int((itera_high - itera_low)/itera_delta) + 1
  itera_size_2 = 1000
  
  allocate(     muu(0:mu_size-1))
  allocate( lyap_mu(0:mu_size-1))
  allocate(   itera(0:itera_size-1))
  allocate(     der(0:itera_size_2-1))
  allocate(      ln(0:itera_size_2-1))
  allocate(  numeri(0:itera_size_2-1))
  allocate(lyapunov(0:itera_size_2-1))
  
  !-----------------------------------------------------------------------------------
  ! Generation of mu coefficients and time values with custom rewrite of
  ! python arange function

  call arange(mu_low, mu_high, mu_delta, muu)
  call arange(itera_low, itera_high, itera_delta, itera)

  x_np1   = 0.0_dp
  lyap_mu = 0.0_dp
  
  open(28, file= "data/lyap_mu.dat")
  
  !-----------------------------------------------------------------------------------
  ! Multiple system evolution for each values of mu,

  mu_cycle:do i = 0, mu_size-1

    !---------------------------------------------------------------------------------
    ! Postion is reinitialized for each mu values

    x_n      = 0.1_dp
    lyapunov = 0.0_dp

    write (file_name,"('data/data_mu_',f19.17,'.dat')") muu(i)
    open(i, file = trim(file_name))
    
    lyapunov_calc:do j = 0, itera_size-1

      numeri    = 0.0_dp
      der       = 0.0_dp
      ln        = 0.0_dp

      numeri(0) = x_n
      der(0)    = abs(df(numeri(0),muu(i)))
      ln(0)     = log(der(0))
      
      
      !--------------------------------------------------------------------------------
      ! Evolution via Euler method of the logistic map

      system_evolution:do k=0, int(itera(j))-1

        x_np1     = f(x_n, muu(i))
        numeri(k) = x_np1 
        x_n       = x_np1
          
      end do system_evolution
      
      !--------------------------------------------------------------------------------
      ! Calculation of logarithm for the derivate of the logistic map, needed for
      ! Lyapunov coefficient calculation

      ln_and_der_calc:do l = 0, int(itera(j))-1
      
        der(l) = abs(df(numeri(l),muu(i)))
        ln(l)  = log(der(l))
      
      end do ln_and_der_calc

      log_sum     = sum(ln)
      lyapunov(j) = log_sum/itera(j)

      write(i,*) itera(j), lyapunov(j)
      
    end do lyapunov_calc

    lyap_mu(i) = lyapunov(itera_size-1)
    
    close(i)
    write(28,*) muu(i), lyap_mu(i) 
    
  end do mu_cycle

  !------------------------------------------------------------------------------------
  ! Closing procedures

  close(28)  
  
  deallocate(muu)
  deallocate(lyap_mu)
  deallocate(itera)
  deallocate(der)
  deallocate(ln)
  deallocate(numeri)
  deallocate(lyapunov)

contains


!--------------------------------------------------------------------------------------
! Logistic map definition

real(dp) function f(x_in,mu_in)
  
  use carepackage
  implicit none

  real(dp), intent(in) :: x_in, mu_in

  f = mu_in*x_in*(1.0_dp-x_in)

end function f

!--------------------------------------------------------------------------------------
! Derivative of the logistic map definition, needed for Lyapunov calculation

real(dp) function df(x_in,mu_in)

  use carepackage
  implicit none

  real(dp), intent(in) :: x_in, mu_in

  df = mu_in*(1.0_dp-(2.0_dp*x_in))

end function df

end program
