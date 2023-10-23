SUBROUTINE compute_TRpr(nloci, X)

    IMPLICIT NONE
    INTEGER :: nloci, ngc
    DOUBLE PRECISION, DIMENSION(3,3,3) :: P
    DOUBLE PRECISION, DIMENSION(3**nloci,3**nloci,3**nloci) :: X
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Xc, Xt
    INTEGER :: n,m,i,j, k, kk, kkk, l, ll, lll, z
    
    P(1, 1, :) = (/  1.0  ,0.0   ,0.0  /) ! AA and AA gives AA
    P(1, 2, :) = (/  0.5  ,0.5   ,0.0  /) ! AA and Aa gives 0.5AA + 0.5Aa
    P(1, 3, :) = (/  0.0  ,1.0   ,0.0  /) ! AA and aa gives Aa
    P(2, 1, :) = (/  0.5  ,0.5   ,0.0  /) ! Aa and AA gives 0.5AA + 0.5Aa
    P(2, 2, :) = (/  0.25 ,0.5   ,0.25 /) ! Aa and Aa gives 0.25AA + 0.5Aa + 0.25aa
    P(2, 3, :) = (/  0.0  ,0.5   ,0.5  /) ! Aa and aa gives 0.5Aa + 0.5aa
    P(3, 1, :) = (/  0.0  ,1.0   ,0.0  /) ! aa and AA gives Aa 
    P(3, 2, :) = (/  0.0  ,0.5   ,0.5  /) ! aa and Aa gives 0.5Aa + 0.5aa
    P(3, 3, :) = (/  0.0  ,0.0   ,1.0  /) ! aa and aa gives aa

    IF(nloci == 1)THEN
        X(1:3,1:3,1:3) = P    
    ELSE
        ALLOCATE(Xc(3,3,3)) 
        Xc = P
        DO n=2,nloci
            ngc = 3**n
            ALLOCATE(Xt(ngc,ngc,ngc)) 

            DO i=1,ngc
                k = i-1
                kk = INT(k/3)
                kkk = k - 3*kk
                DO j=i,ngc
                    l = j-1
                    ll = INT(l/3)
                    lll = l - 3*ll
                    z = 1
                    DO m=1,3**(n-1)
                        Xt(i, j, z:(z+2)) = P(1 + kkk, 1 + lll,:) * Xc(1 + kk, 1 + ll,m)
                        z = z + 3
                    END DO
                END DO
                DO j=1,i-1
                    Xt(i, j,:) = Xt(j, i,:)
                END DO
            END DO
            DEALLOCATE(Xc)
            ALLOCATE(Xc(ngc,ngc,ngc))
            Xc = Xt
            DEALLOCATE(Xt)
        END DO
        X = Xc
        DEALLOCATE(Xc)
    END IF

END SUBROUTINE compute_TRpr

