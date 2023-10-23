SUBROUTINE compute_HWpr(allef, nloci, X)

    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(nloci) :: allef
    DOUBLE PRECISION, DIMENSION(3) :: P, Pc
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Xc, Xt
    DOUBLE PRECISION, DIMENSION(3**nloci) :: X
    INTEGER :: nloci, k
    
    P = (/  (1 - allef(1))**2 , 2*allef(1)*(1-allef(1)) , allef(1)**2  /) 

    IF (nloci == 1) THEN
        X(1:3) = P
    ELSE
        ALLOCATE(Xc(3)) 
        Xc = P
        DO k=2,nloci
            Pc = (/ (1 - allef(k))**2, 2 * allef(k) * (1 - allef(k)), allef(k)**2 /)
            ALLOCATE(Xt(3**k)) 
            Xt = PACK(SPREAD(Xc, dim=2, ncopies=3), mask = .true.) * PACK(SPREAD(Pc, dim=1, ncopies=3**(k-1)), mask = .true.)
            DEALLOCATE(Xc)
            ALLOCATE(Xc(3**k))
            Xc = Xt
            DEALLOCATE(Xt)
        END DO 
        X = Xc
        DEALLOCATE(Xc)
     END IF
 
END SUBROUTINE compute_HWpr
