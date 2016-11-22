module tableio
IMPLICIT NONE

contains

subroutine datainput(fname,table,length)
IMPLICIT NONE
	character(80),intent(in) :: fname
	integer,intent(out) :: length
	real(8),dimension(:,:),allocatable,intent(out) :: table
	
	integer :: i, ios
	
	open(unit = 11, file = fname, status = 'old', action = 'read', iostat = ios)
	if (ios .ne. 0) then
		write(*,'(a,i3,a,a)') 'error opening unit', 11, ' -- ', fname
		stop('END PROGRAM')
	endif
	i = 0
	do
		i = i + 1
		read(11,*,iostat = ios)
		if (ios .lt. 0) then
			length = i - 1
			rewind(11)
			exit
		endif
	enddo
	allocate(table(length,2))
	do i = 1,length
		read(11,*) table(i,:)
	enddo
	close(11)
endsubroutine datainput

endmodule tableio