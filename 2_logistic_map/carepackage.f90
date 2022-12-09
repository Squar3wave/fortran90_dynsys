module carepackage

  integer, parameter, public :: dp = 8

contains  
  
subroutine arange(low_bound, hi_bound, range_arr, array)

  real(dp), intent(in) :: low_bound, hi_bound, range_arr
  real(dp), intent(out) :: array(0:)
  integer :: n, i
    
  !n = size(array)
   n =int((hi_bound - low_bound)/range_arr) + 1 
 
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
  
end module carepackage
