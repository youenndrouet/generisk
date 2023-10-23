SUBROUTINE likproc(pedbrut, n, ng, ndis, nloci, ft, dis, age, tes, asc, liknum, likdenom)

  IMPLICIT NONE
  INTEGER, DIMENSION(n,5 + nloci + 2*ndis), INTENT(IN) :: pedbrut 
  DOUBLE PRECISION, DIMENSION(121,2,ndis,ng), INTENT(IN) :: ft
  INTEGER, DIMENSION(n,ndis), INTENT(IN) :: age, dis
  DOUBLE PRECISION, DIMENSION(n,ndis,ng) :: aff, una
  DOUBLE PRECISION, DIMENSION(n,ng), INTENT(OUT) :: liknum, likdenom
  INTEGER, DIMENSION(n,ng), INTENT(IN) :: tes, asc
  INTEGER, INTENT(IN) :: n, ng, ndis, nloci
  INTEGER :: i,j,c
  
  DO i=1,n !loop over individuals
     DO j=1,ndis !loop over diseases

       IF(pedbrut(i,2) == 0) THEN !female
          aff(i,j,:) = ft(age(i,j)+1, 2,j,:) !+1 because ft begins at age 0, last dim is genotypes
       ELSE !male
          aff(i,j,:) = ft(age(i,j)+1, 1,j,:) 
       END IF

       ! special cases
       c = 5 + nloci + 2*j
       SELECT CASE(pedbrut(i, c-1)) 
         CASE(-1)! UNKNOWN phenotypic status 
           aff(i,j,:) = (/ 1, 1, 1 /) ! aff == 1 (tric to not influence the likelihood) gives same result  than dis(i)=0; aff(i,:) = (/ 0, 0, 0 /) 
         CASE(0) ! unaffected  
           IF(pedbrut(i, c) == 0) THEN ! UNKNOWN age indicated by 0 value
              aff(i,j,:) = (/ 0, 0, 0 /)
           END IF
         CASE(1) ! affected 
           IF(pedbrut(i, c) == 0) THEN ! UNKNOWN age 
              aff(i,j,:) = (/ 1, 1, 1 /)
           END IF
       END SELECT 

     END DO  
  END DO
  
  una = 1 - aff
  
  DO i=1,ng !loop over genotypes 
     !As we assume unrelated diseases, likelihood is simply the PRODUCT over individual diseases risks
     liknum(:,i)    = PRODUCT(dis*aff(:,:,i) + (1 - dis)*una(:,:,i), dim =2) * tes(:,i)      ! Pr(P,G)
     likdenom(:, i) = PRODUCT(dis*aff(:,:,i) + (1 - dis)*una(:,:,i), dim =2) * asc(:,i)      ! Pr(P,Gindex)
  END DO

END SUBROUTINE likproc
