module carepackage

  integer, parameter, public :: dp = 8

contains  

!---------------------------------------------------------------------------------------
! Universal use

subroutine arange(low_bound, hi_bound, range_arr, array)

  real(dp), intent(in) :: low_bound, hi_bound, range_arr
  real(dp), intent(out) :: array(0:)
  integer :: n, i
    
  n = size(array)
 
  do i=0, n-1
      
    if (n == 0) return
      
    if (n == 1) then
      
      array(0) = low_bound
      array(1) = hi_bound
      
      return
      
    end if
      
    if (i==0) then
      
      array(i) = low_bound
      
    elseif (i>0) then 
      
      array(i) = array(i-1) + range_arr
      
    end if
    
  end do
    
end subroutine

subroutine linspace(low_bound, hi_bound, array)

  real(dp), intent(in) :: low_bound, hi_bound
  real(dp), intent(out) :: array(:)
  real(dp) :: range_arr
  integer :: n, i
  n = size(array)
  range_arr = hi_bound - low_bound

  if (n == 0) return

  if (n == 1) then

    array(1) = low_bound
    return

  end if


  do i=1, n
  
    array(i) = low_bound + range_arr * (i - 1) / (n - 1)
  
  end do

end subroutine

subroutine lin_fit(x_in, y_in,size_arr, coeff_out)

  real(dp), intent(in)     :: x_in(0:), y_in(0:)
  integer, intent(in)      :: size_arr
  integer                  :: i,j
  real(dp), intent(out)    :: coeff_out(0:1)
  real(dp), dimension(:)   :: coeff(0:1), x_1_temp(0:size_arr-1), x_2_temp(0:size_arr-1)
  real(dp), dimension(:,:) :: matrix(0:1, 0:1), matrix_inv(0:1, 0:1)
  real(dp)                 :: invert
  
  do i=0,1
   do j=0,1
   
     x_1_temp    = x_in**i
     x_2_temp    = x_in**j
   
     matrix(i,j) = dot_product(x_1_temp,x_2_temp)
   
   end do
  end do

  invert = 1.0_dp/(matrix(0,0)*matrix(1,1) - matrix(0,1)*matrix(1,0))
  
  matrix_inv(0,0) =  matrix(1,1)
  matrix_inv(0,1) = -matrix(0,1)   
  matrix_inv(1,0) = -matrix(1,0)
  matrix_inv(1,1) =  matrix(0,0)
   
  matrix_inv      = invert * matrix_inv
   
   do i=0,1

     x_1_temp = x_in**i
  
     coeff(i) = dot_product(x_1_temp, y_in)
   
   end do
   
   coeff_out = matmul(matrix_inv,coeff)
   
end subroutine
  
subroutine array_max_finder(array,max_found)

  real(dp), dimension(:), intent(in) :: array(0:)
  real(dp), intent(out) :: max_found
  integer :: i, arrsize

  arrsize=size(array)
  max_found = array(0)

  do i=1, arrsize-1

    if(max_found < array(i)) then

      max_found = array(i)

    else

      max_found = max_found
  
    end if

  end do

end subroutine

!---------------------------------------------------------------------------------------
! Specific use (define here if needed)
  
end module carepackage
