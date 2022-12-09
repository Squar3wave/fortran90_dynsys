program albero

  use carepackage

  IMPLICIT NONE

  integer  :: n_scala, n_tot, i, j, k, q, nn, nk
  real(dp) :: delta, gg, phi, a_i
  real(dp), allocatable, dimension (:)   :: p, scales, ln_scales
  real(dp), allocatable, dimension (:,:) :: moment, ln_mom
  real(dp), dimension (:,:)              :: lin_coeff(0:1,0:3)                                            

  open( 0, file="signal.dat")
  open( 1, file="mft_mom_1.dat")
  open( 2, file="mft_mom_2.dat")
  open( 3, file="mft_mom_3.dat")  
  open(11, file="mft_lnm_1.dat")
  open(12, file="mft_lnm_2.dat")
  open(13, file="mft_lnm_3.dat")
  
  n_scala = 15
  delta   = 0.8_dp
  n_tot   = 2**n_scala

  allocate(        p(0:n_tot-1))
  allocate(   scales(0:n_scala-1))
  allocate(ln_scales(0:n_scala-1))
  allocate(   moment(0:n_scala-1,0:3))
  allocate(   ln_mom(0:n_scala-1,0:3))

  scales = 0.0_dp
  moment = 0.0_dp
  
  
  !------------------------------------------------------------------
  ! first we set all values to 1

  p = 1.0_dp

  !------------------------------------------------------------------
  ! Ladder and signal generation

  do i = 0, n_scala-1

    !----------------------------------------------------------------  
    ! Ladder parameter definition
  
    nn        = 2**i
    nk        = int(n_tot/nn)
    scales(i) = (2.0_dp)**(-i)
    
    !----------------------------------------------------------------
    ! Generation of the signal, which biforks more and more
    ! as i grows in value
    
    do j = 0, nn-1
  
      call random_number(gg)
  
      phi = 1.0_dp - delta + 2.0_dp*delta*gg
    
      do k = j*nk, (j+1)*nk-1
    
        p(k) = p(k)*phi
    
      end do
  
    end do
    
    !----------------------------------------------------------------
    ! Since we're dealing with a power law we calculate all 
    ! needed natural logarithms, in order to later plot a more
    ! approachable linear fit 
    
    ln_scales(i) = log(scales(i))
  
  end do

  !------------------------------------------------------------------  
  ! Writing the results on file on a separate loop in order to
  ! only store the final ones
  
  do k = 0, size(p)-1
    
    write(0,*) k, p(k)
    
  end do
  
  !------------------------------------------------------------------  
  ! We deal with the moments calculation in a separate do loop 
  ! in order to avoid using temporary values
  
  do i = 0, n_scala-1

    !----------------------------------------------------------------  
    ! Ladder parameter definition
  
    nn = 2**i
    nk = n_tot/nn

    !----------------------------------------------------------------
    ! 
     
    do j = 0, nn-1
  
      a_i = 0
    
      do k = j*nk, (j+1)*nk-1
    
        a_i = a_i + p(k)/nk
    
      end do
  
      do q = 0, 3
    
        moment(i,q) = moment(i,q) + (a_i**q)/nn
    
      end do
      
    end do

  end do

  do i=0,n_scala-1
  
    !-------------------------------------------------------------
    ! Since we're dealing with a power law we calculate all 
    ! needed natural logarithms, in order to later plot a more
    ! approachable linear fit 
        
    ln_mom(i,1) = log(moment(i,1))
    ln_mom(i,2) = log(moment(i,2))
    ln_mom(i,3) = log(moment(i,3))    
    
  end do
  
  !------------------------------------------------------------------  
  ! 
  
  call lin_fit(ln_scales,ln_mom(:,1),n_scala,lin_coeff(:,1))
  call lin_fit(ln_scales,ln_mom(:,2),n_scala,lin_coeff(:,2))
  call lin_fit(ln_scales,ln_mom(:,3),n_scala,lin_coeff(:,3)) 

  print "('-----------------------------------------------')"
  print "('| Linear fit results for the three momentums  |')"
  print "('|---------------------------------------------|')"
  print "('|       |    offset    |  angular coefficient |')"
  print "('|-------|--------------|----------------------|')"
  print "('| mom_1 | ', f11.7, '  |  ', f11.7, '         |')", lin_coeff(:,1)
  print "('| mom_2 | ', f11.7, '  |  ', f11.7, '         |')", lin_coeff(:,2)
  print "('| mom_3 | ', f11.7, '  |  ', f11.7, '         |')", lin_coeff(:,3)
  print "('-----------------------------------------------')"

  do i=0,n_scala-1
  
    write( 1,*)    scales(i), moment(i,1)
    write( 2,*)    scales(i), moment(i,2)
    write( 3,*)    scales(i), moment(i,3)  
    write(11,*) ln_scales(i), ln_mom(i,1)
    write(12,*) ln_scales(i), ln_mom(i,2)
    write(13,*) ln_scales(i), ln_mom(i,3)  
    
    end do
  
  close( 0)
  close( 1)
  close( 2)
  close( 3)
  close(11)
  close(12)
  close(13)
   
end program
