PROGRAM conversion

  IMPLICIT NONE
 
  INTEGER :: n1, n2, i, j
  CHARACTER(LEN=250) :: arg1, arg2, arg3, arg4
  CHARACTER(LEN=500) :: file1, file2
  CHARACTER(LEN=*), PARAMETER :: s1="file1", s2="file2"
  REAL:: x, y
  REAL, DIMENSION(:,:), ALLOCATABLE :: zz

  CALL GETARG(1,arg1)
  CALL GETARG(2,arg2)
  CALL GETARG(3,arg3) 
  CALL GETARG(4,arg4) 

  READ(arg1,*) n1
  READ(arg2,*) n2

  ALLOCATE(zz(n1,n2))

  file1 = TRIM(ADJUSTL(arg3))//s1//TRIM(ADJUSTL(arg4))
  file2 = TRIM(ADJUSTL(arg3))//s2//TRIM(ADJUSTL(arg4))


  OPEN (1,file=file1, status="old")

  DO i=1,n1
     DO j=1,n2
        read(1,*,END=10) x,y,zz(i,j)
     END DO
  END DO
10 CLOSE(1)


  OPEN (2,file=file2,form='unformatted')
  WRITE(2)n1,n2
  WRITE(2)zz
  CLOSE(2)

  
END PROGRAM conversion
