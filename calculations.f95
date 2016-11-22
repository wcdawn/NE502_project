module calculations
IMPLICIT NONE

contains
subroutine interpolate(table,length,x,y)
	IMPLICIT NONE
	real(8),dimension(:,:),intent(in) :: table
	real(8),intent(in) :: x
	real(8),intent(out) :: y
	
	integer,intent(in) :: length
	integer :: i
	
	integer :: i_hi, i_lo
	real(8) :: r_hi, r_lo
	real(8) :: slope
	
	101 format(a)
	103 format(a,e12.6)
	
	if ((x .gt. table(length,1)) .or. (x .lt. table(1,1))) then
		write(*,101) 'x value out of bounds for interpolation'
		write(*,103) 'x = ', x
		write(*,103) 'x_max = ', table(length,1)
		write(*,103) 'x_min = ', table(1,1)
		STOP
	endif
	
	! table in order from lowest x to highest x
	! x is in first column of 2 column table
	do i = 1,length
		if (table(i,1) .gt. x) then
			i_hi = i
			i_lo = i - 1
			r_hi = table(i_hi,1)
			r_lo = table(i_lo,1)
			exit
		endif
	enddo
	
	slope = (table(i_hi,2) - table(i_lo,2)) / (table(i_hi,1) - table(i_lo,1))
	y = slope * (x - table(i_hi,1)) + table(i_hi,2)
endsubroutine interpolate

! real(8) function qpp(z)
subroutine forwarddiff(array_in,length,delta_z,array_out)
	IMPLICIT NONE
	integer,intent(in) :: length
	real(8),dimension(:),intent(in) :: array_in
	real(8),intent(in) :: delta_z
	real(8),dimension(:),allocatable,intent(out) :: array_out
	
	integer :: i
	
	allocate(array_out(length - 1))
	do i = 1,length-1
		array_out(i) = (array_in(i + 1) - array_in(i)) / delta_z
	enddo
	
endsubroutine forwarddiff

subroutine trapz(array_in,length,mesh,low,high,sol)
	IMPLICIT NONE
	integer,intent(in) :: length
	real(8),dimension(:),intent(in) :: array_in, mesh
	real(8),intent(in) :: low, high
	real(8),intent(out) :: sol
	
	integer :: i
	integer :: lo_lo, lo_hi, hi_lo, hi_hi
	real(8) :: array_low, array_high

	101 format(a)
	103 format(a,e12.6)
	
	if (((low .lt. mesh(1))) .or. (high .gt. mesh(length))) then
		write(*,101) 'value out of bounds for integration'
		STOP
	endif	
	
	! assume mesh from low to high
	do i = 1,length
		if (mesh(i) .gt. low) then
			lo_hi = i
			lo_lo = i - 1
			exit
		endif
	enddo
	do i = 1,length
		if (mesh(i) .gt. high) then
			hi_hi = i
			hi_lo = i - 1
			exit
		endif
	enddo
	
	
	! first and last trap will be different		
	sol = 0.0d0
	! first
	array_low = ((array_in(lo_hi) - array_in(lo_lo)) / (mesh(lo_hi) - mesh(lo_lo))) * (low - mesh(lo_lo)) + array_in(lo_lo)
	sol = sol + (5.0d-1) * (mesh(lo_hi) - low) * (array_in(lo_hi) + array_low)
		! general
	do i = lo_hi,hi_lo
		sol = sol + (5.0d-1) * (mesh(i + 1) - mesh(i)) * (array_in(i + 1) + array_in(i))
	enddo
	! last
	array_high = ((array_in(hi_hi) - array_in(hi_lo)) / (mesh(hi_hi) - mesh(hi_lo))) * (high - mesh(hi_lo)) + array_in(hi_lo)
	sol = sol + (5.0d-1) * (mesh(hi_hi) - high) * (array_in(hi_hi) + array_high)
	
endsubroutine trapz





endmodule calculations