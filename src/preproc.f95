
SUBROUTINE preproc(n,ng,nloci,ndis,pedbrut,ped,kmax,ids,dis,tes,asc,age)
     	 	 
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: n      !pedigree size
     INTEGER, INTENT(IN) :: ng     !number of genotypes (3**nloci)
     INTEGER, INTENT(IN) :: nloci  !number of loci
     INTEGER, INTENT(IN) :: ndis   !number of diseases
     INTEGER, INTENT(IN) :: kmax   !max number of child/spouse/siblings for one individual
     INTEGER, DIMENSION(n,5 + nloci + 2*ndis), INTENT(IN) :: pedbrut

     INTEGER, DIMENSION(n,5,kmax), INTENT(OUT) :: ids
     INTEGER, DIMENSION(n,4), INTENT(OUT) :: ped
     INTEGER, DIMENSION(n,ng), INTENT(OUT) :: tes, asc
     INTEGER, DIMENSION(n,ndis), INTENT(OUT) :: dis, age

     INTEGER, DIMENSION(n) :: newid, newfather, newmother
     INTEGER :: i
	 
     ids    = 0  ! table of individual id specifying all genetic links between relatives (keeps in memory the identifiers of siblings, spouses and childs of all individuals)

     ped(:,1:4) = pedbrut(:,1:4) !--pedbrut must be an array with first columns 1:id | 2:gender | 3:fatherid | 4:motherid
     !--recoding of ids to simplify implemention
     newid     = (/(i, i=1, n)/)
     newfather = 0
     newmother = 0

    DO i=1,n

        IF (ped(i,3) /= 0) THEN
            newfather(i) = posint(ped(i,3), ped(:,1), n)
        ELSE
            newfather(i) = 0
        END IF
        IF (ped(i,4) /= 0) THEN
            newmother(i) = posint(ped(i,4), ped(:,1), n)
        ELSE
            newmother(i) = 0
        END IF

    END DO

    ped(:,1) = newid
    ped(:,3) = newfather
    ped(:,4) = newmother
    ! now individuals id indicate the respectives position of the individuals in the pedigree table

    ! code above find childs, spouses and siblings of all individuals and keep this in memory for future implementation
	CALL findids
	! code above prepare objects for likelihood calculation
	CALL preplik
	 

  CONTAINS 
	 
	 !--- findids
		 
	SUBROUTINE findids() 

	  IMPLICIT NONE
	  
	  INTEGER :: i, j, mid, fid
	  INTEGER, DIMENSION(kmax) :: cid, sid, fsibid, sfid, smid 
	  INTEGER :: kc, kfs, ksf, ksm 
		
		! --cid : childs ids 
		! --sid : spouses ids
		! --fsibid : siblings ids
		! --sfid : mother's ids of half sibling (ie sib from father)
		! --smid : father's ids of half sibling (ie sib from mother)
		
		DO i=1,n
		
		  cid = 0
		  sid = 0
		  fsibid = 0
		  sfid = 0
		  smid = 0 
		  fid = ped(i,3)
		  mid = ped(i,4)
		  
		  kc = 0; kfs = 0; ksf = 0; ksm = 0 
		  
		  IF ((mid /= 0).OR.(fid /= 0)) THEN 
			
			! search siblings
			IF ((mid /= 0).AND.(fid /= 0)) THEN
				DO j=1,n
				  IF((ped(j,4) == mid).AND.(ped(j,3) /= fid)) THEN   
					 ksm = ksm + 1
					 smid(ksm) = ped(j,3)
				  ELSE IF ((ped(j,4) /= mid).AND.(ped(j,3) == fid)) THEN  
					 ksf = ksf + 1 
					 sfid(ksf) = ped(j,4)
				  END IF
				END DO
			ELSE 
				IF (mid /= 0) THEN
				  DO j=1,n
					IF((ped(j,4) == mid).AND.(ped(j,3) /= fid)) THEN   
					  ksm = ksm + 1
					  smid(ksm) = ped(j,3)
					ELSE IF (ped(j,4) == -1) THEN  
					  ksf = ksf + 1 
					  sfid(ksf) = ped(j, 3)
					END IF
				  END DO
				ELSE 
				  DO j=1,n
					IF (ped(j,4) == -1) THEN   
					  ksm = ksm + 1
					  smid(ksm) = ped(j,3)
					ELSE IF ((ped(j,3) == fid).AND.(ped(j,4) /= mid)) THEN
					  ksf = ksf + 1 
					  sfid(ksf) = ped(j,4)
					END IF
				  END DO
				END IF
			END IF
		  
			!search full-siblings
			DO j=1,n 
			  IF((ped(j,3) == fid).AND.(ped(j,4) == mid).AND.(ped(j,1) /= i)) THEN
				kfs = kfs + 1
				fsibid(kfs) = ped(j, 1)
			  END IF
			END DO
			
		  ELSE 
			
			DO j=1,n 
			  IF(ped(j,4) == -1) THEN
				ksm = ksm + 1
				smid(ksm) = ped(j,3)
			  END IF
			  IF(ped(j,3) == -1) THEN
				ksf = ksf + 1 
				sfid(ksf) = ped(j,4)
			  END IF
			  IF((ped(j,3) == -1).AND.(ped(j,4) == -1))THEN
				kfs = kfs + 1
				fsibid(kfs) = ped(j,1)
			  END IF
			END DO
			  
		  END IF
		  
		 ! loop to find children and spouses
		  IF (ped(i, 2) == 1) THEN ! indid i is a male, so we look for a wife
			  DO j=1,n
				IF (ped(j,3) == i) THEN !found a child
				   kc = kc + 1
				   cid(kc) = j
				   sid(kc) = ped(j,4)
				END IF
			  END DO
		  ELSE                    ! indid i is a female, so we look for a husband
			  DO j=1,n
				IF (ped(j,4) == i) THEN !found a child
				   kc = kc + 1
				   cid(kc) = j
				   sid(kc) = ped(j,3)
				END IF
			  END DO
		  END IF

		  ids(i,1,:) = smid
		  ids(i,2,:) = sfid
		  ids(i,3,:) = fsibid
		  ids(i,4,:) = sid
		  ids(i,5,:) = cid
		  
		END DO

	END SUBROUTINE findids


	SUBROUTINE preplik()

	  IMPLICIT NONE
	  INTEGER, DIMENSION(:), ALLOCATABLE ::  tescm, tesct
	  INTEGER, DIMENSION(3) :: tesc1
	  
	  INTEGER :: i,j,k,c
	   
	  asc = 1
	  tes = 1
	  
	  DO i=1,n !loop over individuals
	 
		 ! dis (disease indicator) and age (at diagnosis or censory) implementation
		  DO j=1,ndis !loop over diseases
		   c = 5 + nloci + 2*j
		   SELECT CASE(pedbrut(i, c-1)) 
			 CASE(-1)! UNKNOWN phenotypic status 
			   dis(i,j) = 1 
			 CASE(0) ! unaffected  
			   dis(i,j) = 0
			 CASE(1) ! affected 
			   dis(i,j) = 1
		   END SELECT 
		   age(i,j) = pedbrut(i, c)
		  END DO  

		 ! tes (genetic tests results) implementation
		
		 !first loci
		   c = 5+1
		   SELECT CASE(pedbrut(i, c)) 
			 CASE(0)
			  tesc1 = (/ 1, 0, 0 /) !aa
			 CASE(1)
			  tesc1 = (/ 0, 1, 0 /) !Aa
			 CASE(2)
			  tesc1 = (/ 0, 0, 1 /) !AA
			 CASE(4)
			  tesc1 = (/ 1, 1, 1 /) !UNKNOWN, default value
		   END SELECT
		 
		   IF (nloci == 1) THEN
			tes(i,:) = tesc1
		   ELSE

			 !other loci
			 ALLOCATE(tescm(3))
			 tescm = tesc1

			 DO k=2,nloci !loop over loci

				c = 5 + k
				SELECT CASE(pedbrut(i, c)) 
				  CASE(0)
					tesc1 = (/ 1, 0, 0 /) !aa
				  CASE(1)
					tesc1 = (/ 0, 1, 0 /) !Aa
				  CASE(2)
					tesc1 = (/ 0, 0, 1 /) !AA
				  CASE(4)
					tesc1 = (/ 1, 1, 1 /) !UNKNOWN, default value
				END SELECT
							
				ALLOCATE(tesct(3**k)) 
		  tesct = PACK(SPREAD(tescm, dim=2, ncopies=3), mask = .true.) * PACK(SPREAD(tesc1, dim=1, ncopies=3**(k-1)), mask = .true.)
				DEALLOCATE(tescm)
				ALLOCATE(tescm(3**k))
				tescm = tesct
				DEALLOCATE(tesct)

			 END DO
			 tes(i,:) = tescm
			 DEALLOCATE(tescm)

		 END IF

		 IF(pedbrut(i, 5) == 1) THEN !index case
			 asc(i,:) = tes(i,:)
		 END IF
	  END DO
	  

	END SUBROUTINE preplik
	 
	 
	 !--- posint

	 INTEGER FUNCTION posint(x,table,n)
	    !IMPLICIT NONE
		! returns the position of the (first) element of table that is equal to x 
		INTEGER, DIMENSION(n), INTENT(IN)  :: table
		INTEGER, INTENT(IN) :: x, n
		LOGICAL :: fin 
		posint = 1
		fin = .FALSE.
		DO WHILE ((.NOT. fin) .AND. (posint <= n))
		  if (table(posint)==x) THEN
			fin = .TRUE.
		  ELSE
			posint = posint + 1
		  END IF
		END DO 
		IF (.NOT. fin) THEN
		   posint = 0
		END IF
	END FUNCTION posint
	 

END SUBROUTINE preproc



