! FORTRAN code to peel a pedigree by ES algorithm

  SUBROUTINE peeling(loglik, n, ng, kmax, ped, ids, lik, genepr, weight, counsid, hwpr, gpr)

     IMPLICIT NONE

     DOUBLE PRECISION, INTENT(OUT) :: weight, loglik
     DOUBLE PRECISION, DIMENSION(ng), INTENT(OUT) ::genepr ! posterior genotype proba for counselee

     INTEGER, INTENT(IN) :: n      !pedigree size
     INTEGER, INTENT(IN) :: ng
     INTEGER, INTENT(IN) :: kmax
     INTEGER, DIMENSION(n,5,kmax), INTENT(IN) :: ids
     INTEGER, DIMENSION(n,4), INTENT(IN) :: ped
     DOUBLE PRECISION, DIMENSION(n,ng), INTENT(IN) :: lik
     DOUBLE PRECISION, DIMENSION(ng), INTENT(IN) :: hwpr
     DOUBLE PRECISION, DIMENSION(ng,ng,ng), INTENT(IN) :: gpr
     INTEGER, INTENT(IN) :: counsid

     LOGICAL, PARAMETER :: verbose = .FALSE.
     DOUBLE PRECISION, DIMENSION(kmax,ng) :: pp
     CHARACTER(LEN=10), DIMENSION(kmax) :: ppid
     INTEGER :: ppk

     ! init of Global variables
     weight = 0  !weight (offset) used for the pedigree likelihood calculation
     ppk    = 0  !number of posterior() calls
     ppid   = '' !identifiers of posterior() calls
     pp     = 1  !stockage matrice of posterior probabilities (outputs of posterior() calls)

     !--ped must be an array with columns 1:id | 2:gender | 3:fatherid | 4:motherid
     !--above code assume ped is sorted by id, with id=1:n

     CALL peel ! modifies genepr and weight
     loglik = LOG(SUM(genepr)) + weight
     genepr = genepr/SUM(genepr)
	 
	 CONTAINS 
	 
!--- peel
		  
		SUBROUTINE peel() 
		! FORTRAN code for the Elston-Stewart peeling function
		! which calculate the likelihood of pedigrees without loop
		! -- REF : Fernando RL, Stricker C, Elston RC. An efficient 
		!algorithm to compute the posterior genotypic distribution 
		!for every member of a pedigree without loops. Theor. Appl. Genet. 1993, 87:89-93.

		! largely inspired from the peeling function of the BayesMendel R package

		 IMPLICIT NONE
		 
		 DOUBLE PRECISION, DIMENSION(ng) :: a, p, temp
		 INTEGER :: i,k
		 INTEGER, DIMENSION(kmax) :: sid
		   
		 IF (verbose) THEN 
			OPEN(UNIT=1, FILE="debug.log", FORM = 'formatted')
			CALL debugfile
		 END IF
		 
		  a = 1;  p = 1; 
				
			sid = ids(counsid,4,:)
			CALL uid(sid,kmax) ! sid has unique values  
			IF (verbose) THEN 
			   WRITE(1,'(A)') 'tracer : peel() calling posterior()'
			END IF

			IF (SUM(sid) > 0)  THEN ! counselee had at least a child
				k = 1
				DO WHILE(sid(k) /= 0) ! loop over the children
					IF (verbose) THEN 
					  WRITE(1,'(A)') 'tracer : peel() calling posterior()'
					END IF
					CALL posterior(counsid, sid(k)) !modify pp (table of posterior prob)
					p = p * (pp(ppk,:)/SUM(pp(ppk,:)))
					weight = weight + LOG(SUM(pp(ppk,:)))
					k = k + 1       
					IF (verbose) THEN 
						WRITE(1,'(A, 3ES11.3)') 'p', p
					END IF
				END DO
			END IF
			
			IF (verbose) THEN 
			  WRITE(1,'(A)') 'tracer : peel() calling anterior()'
			END IF
			CALL anterior(a, counsid) ! modify a
			IF (verbose) THEN 
			   WRITE(1,'(A, 3ES11.3)') 'p', p
			END IF
			genepr = a/SUM(a) * lik(counsid,:)/SUM(lik(counsid,:)) * p/SUM(p)
			weight = weight + LOG(SUM(a)) + LOG(SUM(lik(counsid,:))) + LOG(SUM(p)) 

			IF (verbose) THEN 
				!write the likelihood results 
				IF(ppk /= 0) THEN
				  WRITE(1,'(A)') 'Posterior (pp) = '
				  DO i=1,ppk
					 WRITE(1,'(A)') ppid(i)
					  WRITE(1,'(3ES11.3)') pp(i,:) 
				  END DO
				END IF
		  
				WRITE(1,'(A)') 'weight = '
				WRITE(1,'(1ES11.3)') weight
				CLOSE(UNIT=1)
			END IF
			
		END SUBROUTINE peel

	 
!------ anterior
	 
	 SUBROUTINE anterior(output, indid) 

	  IMPLICIT NONE
	 
	  INTEGER :: j,k,u,ui,uj,i
	  INTEGER, INTENT(IN) :: indid
	  INTEGER, DIMENSION(kmax) :: fsibid
	  DOUBLE PRECISION, DIMENSION(ng, ng) :: fullsibs, fullsibsi, tgfs
	  DOUBLE PRECISION, DIMENSION(ng) :: output, vec_fprtgfs, mothers, fathers, fsibi  
	  
	  INTEGER, DIMENSION(kmax) :: siblings ! ids of siblings of indid
		
		siblings  = 0
		fullsibs  = 1
		fullsibsi = 1
			
		!ped mut be an array with columns 1:id | 2:gender | 3:fatherid | 4:motherid
		IF ((ped(indid, 3) == 0 ).AND.(ped(indid, 4) == 0))THEN
			output = hwpr
			!WRITE(1,'(A, 3ES11.3)') 'ant(from anterior)',ant
		ELSE 
			IF (verbose) THEN 
			  WRITE(1,'(A)') 'tracer : anterior() calling antmfs(mode = "m")'
			END IF

			CALL antmfs(mothers, indid, 'm', 0) !modify mothers 
			!WRITE(1,'(A, 3ES11.3)') 'mothers(from anterior, mod by antmfs)', mothers
			IF (verbose) THEN 
			  WRITE(1,'(A)') 'tracer : anterior() calling antmfs(mode = "f")'
			END IF

			CALL antmfs(fathers, indid, 'f', 0) !modify fathers
			!WRITE(1,'(A, 3ES11.3)') 'fathers(from anterior, mod by antmfs)', fathers

			fsibid = ids(indid,3,:)
			CALL uid(fsibid, kmax) ! fsibid contains the unique values of ids(indid,3,:)
			
			IF (SUM(fsibid) > 0) THEN !at least one sibling
				k = 1

				DO WHILE(fsibid(k) /= 0) !loop over siblings
					IF (verbose) THEN 
					  WRITE(1,'(A)') 'tracer : anterior() calling sibchild(mode="s")'
					END IF

					CALL sibchild(fsibi,indid, fsibid(k), 0, 's') !modify fsibi
					!WRITE(1,'(A, 3ES11.3)') 'fsibi(from anterior, mod by sibchild)', fsibi

					DO ui=1,ng
					  DO uj=ui,ng
						fullsibsi(ui, uj) = SUM(gpr(ui,uj,:) * fsibi)
					  END DO
					END DO

					DO ui=2,ng 
					  DO uj=1,(ui - 1) 
						fullsibsi(ui, uj) = fullsibsi(uj, ui)
					  END DO
					END DO

					fullsibs = fullsibs * fullsibsi
					k = k + 1

				END DO
			END IF
			
			DO u=1,ng 
				!WRITE(1,'(A, 3ES11.3)') 'fullsibs', fullsibs(u,:)
				tgfs = TRANSPOSE(gpr(:,:,u)*fullsibs)
				DO i=1, ng
				   vec_fprtgfs(i) =  DOT_PRODUCT(fathers,tgfs(:,i)) 
				   !WRITE(1,'(A, 3ES11.3)') 'vec_fprtgfs', vec_fprtgfs(i)
				END DO
				output(u) = DOT_PRODUCT(mothers, vec_fprtgfs)
			END DO
			
		END IF

	END SUBROUTINE anterior
	
	
 !------ antmfs
	
    SUBROUTINE antmfs(output, indid, mode, sid) 

	  ! note : sid argument only /= 0 if mode =='s'

	  IMPLICIT NONE

	  CHARACTER(LEN=1) :: mode 
	  DOUBLE PRECISION, DIMENSION(ng) :: am, likm, af, likf, as, liks, postm, postf, posts, p, output
	  INTEGER, DIMENSION(kmax) :: ssid, sfid, smid  
	  INTEGER :: posi, k, fid, mid, sid
	  INTEGER, INTENT(IN) :: indid
		
		am = 1;    af = 1;    as = 1
		likm = 1;    likf = 1;    liks = 1
		postm = 1;    postf = 1;    posts = 1
		
		SELECT CASE(mode)
		
		CASE('m') !mothers
		
			mid = ped(indid,4)
			
			IF (mid /= 0) THEN
				smid = ids(indid,1,:)
				CALL uid(smid, kmax) !smid contains unique sfid
			   IF(SUM(smid) > 0) THEN
				  k = 1
				  DO WHILE(smid(k) /= 0)
					  IF (verbose) THEN 
						WRITE(1,'(A)') 'tracer : antmfs(mode="m") calling posterior(mid,smid(k))'
					  END IF
					  CALL posterior(mid, smid(k)) !modify pp, ppid and ppk
					  p = pp(ppk,:)
					  postm = postm * (p/SUM(p))
					  weight = weight + LOG(SUM(p))
					  k = k + 1
				  END DO
			   END IF
			   
			   IF (verbose) THEN 
				  WRITE(1,'(A)') 'tracer : antmfs(mode="m") calling anterior(am,mid)' 
			   END IF
			   CALL anterior(am, mid) !modify am
			
			   likm = lik(mid,:)
			   output = (am/SUM(am)) * (likm/SUM(likm)) * (postm/SUM(postm))
			   weight = weight + LOG(SUM(am)) + lOG(SUM(likm)) + LOG(SUM(postm)) 
			ELSE
			   output = hwpr
			END IF

			
		CASE('f') !fathers
		
			fid = ped(indid,3)
			IF (fid /= 0) THEN
				sfid = ids(indid,2,:)
				CALL uid(sfid, kmax) !sfid contains unique sfid
				IF(SUM(sfid) > 0) THEN
				  k = 1
				  DO WHILE(sfid(k) /= 0)
					  IF (verbose) THEN
						WRITE(1,'(A)') 'tracer : antmfs(mode="f") calling posterior(fid, sfid(k))'
					  END IF
					  CALL posterior(fid, sfid(k)) !modify pp, ppid and ppk
					  p = pp(ppk,:)
					  postf = postf * (p/SUM(p))
					  weight = weight + LOG(SUM(p))
					  k = k + 1
				  END DO
				END IF
				
				IF (verbose) THEN
					WRITE(1,'(A)') 'tracer : antmfs(mode="f") calling anterior(af,fid)'
				END IF
				CALL anterior(af, fid) !modify af
			 
				likf = lik(fid,:)
				output = (af/SUM(af)) * (likf/SUM(likf)) * (postf/SUM(postf))
				weight = weight + LOG(SUM(af)) + lOG(SUM(likf)) + LOG(SUM(postf))    
			ELSE
			   output = hwpr
			END IF
		
		CASE('s') !spouses
		
			ssid = ids(sid,4,:)

			CALL uid(ssid, kmax) !ssid contains unique sid
			
			IF(inint(indid, ssid, kmax)) THEN
				posi = posint(indid, ssid, kmax)
				IF(posi /= kmax) THEN
				   ssid(posi:(kmax - 1)) = ssid((posi + 1):kmax)
				   ssid(kmax) = 0
				END IF
			END IF !ssid does not contain indid 
			
			IF(SUM(ssid) > 0) THEN
				k = 1
				DO WHILE(ssid(k) /= 0)
				   IF (verbose) THEN
					  WRITE(1,'(A)') 'tracer : antmfs(mode="s") calling posterior(sid, ssid(k))' 
				   END IF
				   CALL posterior(sid, ssid(k)) !modify pp, ppid and ppk
					p = pp(ppk,:)
					posts = posts * (p/SUM(p))
					weight = weight + LOG(SUM(p))
					k = k + 1
				END DO
			END IF
			
			IF (verbose) THEN
			  WRITE(1,'(A)') 'tracer : antmfs(mode="s") calling anterior(as,sid)'
			END IF
			CALL anterior(as, sid) !modify as
			liks = lik(sid,:)
			output = (as/SUM(as)) * (liks/SUM(liks)) * (posts/SUM(posts))
			weight = weight + LOG(SUM(as)) + lOG(SUM(liks)) + LOG(SUM(posts))
		END SELECT
		
    END SUBROUTINE antmfs

!--- posterior

	SUBROUTINE posterior(indid, sid) 

	  IMPLICIT NONE
	 
	  INTEGER :: indid, sid, i,k,l,ui,uj
	  DOUBLE PRECISION, DIMENSION(ng) :: pijui, spouses, ci
	  DOUBLE PRECISION, DIMENSION(ng,ng) :: childi, children
	  INTEGER, DIMENSION(kmax) :: allsid, allcid, cid
	  
	  LOGICAL :: inpp
	  
		childi = 1 ;    children = 1; 
		! this is a call to posterior function, so we indent ppk
		ppk = ppk + 1
		WRITE(ppid(ppk),'(I3,A,I3)') indid,' ',sid 
	   
		IF (ppk > 1) THEN
		  inpp = .FALSE. 
		  i = 1
		  DO WHILE ((.NOT.(inpp)).AND.(i < ppk))
			inpp = LGE(ppid(ppk),ppid(i))
			IF(inpp)THEN !char strings are identical
			  ! should not append !
			  !rtn = pp(i,:)
			END IF
			i = i + 1
		  END DO      
		END IF
		
		IF (sid /= 0) THEN
			IF (verbose) THEN 
			  WRITE(1,'(A)') 'tracer : posterior() calling antfms(mode="s")'  
			END IF
			CALL antmfs(spouses, indid, 's', sid) !modify spouses
		ELSE
			spouses = hwpr 
		END IF
			   
		allsid = ids(indid,4,:)
		allcid = ids(indid,5,:)
		cid = 0;    k = 1;    l = 1
		
		DO WHILE (allsid(k) /= 0)
		  IF(allsid(k) == sid) THEN
			cid(l) = allcid(k)
			l = l + 1
		  END IF
		  k = k + 1
		END DO

		IF (SUM(cid) > 0) THEN
			k = 1
			DO WHILE (cid(k) /= 0)
				IF (verbose) THEN 
				  WRITE(1,'(A)') 'tracer : posterior() calling sibchild(mode="c")'        
				END IF
				CALL sibchild(ci, indid, 0, cid(k), 'c') !modify ci
				!WRITE(1,'(A, 3ES11.3)') 'ci', ci

				DO ui=1,ng
					DO uj=1,ng
					  childi(ui, uj) = SUM(gpr(ui,uj,:) * ci)
					END DO
				END DO
				
				children = children * childi
				k = k + 1
				
			END DO
		END IF
		
		IF (ped(indid, 2) == 1) THEN
			pijui = MATMUL(spouses, children)
		ELSE
			pijui = MATMUL(spouses, TRANSPOSE(children))
		END IF
		
		!WRITE(1,'(A, 3ES11.3)') 'pijui', pijui
		pp(ppk,:) = pijui
		
	END SUBROUTINE posterior

!--- sibchild

	SUBROUTINE sibchild(output, indid, fsibid, cid, mode) 
	  
	  IMPLICIT NONE

	  CHARACTER(LEN=1) :: mode 
	  INTEGER, DIMENSION(kmax) :: ssibj, scidj
	  INTEGER :: k, indid, fsibid, cid
	  
	  DOUBLE PRECISION, DIMENSION(ng) :: posts, p, likk, output
	  
		posts = 1
		
		SELECT CASE(mode)
		
		CASE('s')
		
			ssibj = ids(fsibid,4,:)
			CALL uid(ssibj, kmax)         

			IF (SUM(ssibj) > 0) THEN
			  k = 1
			  DO WHILE(ssibj(k) /= 0)
				  IF (verbose) THEN 
					WRITE(1,'(A)') 'tracer : sibchild(case="s") calling posterior()'
				  END IF
				  CALL posterior(fsibid, ssibj(k)) !modify pp, ppid and ppk
				  p = pp(ppk,:)
				  posts = posts * (p/SUM(p))
				  weight = weight + LOG(SUM(p))
				  k = k + 1
			  END DO
			END IF
			likk = lik(fsibid,:)
			output = (likk/SUM(likk)) * (posts/SUM(posts))
			weight = weight + LOG(SUM(likk)) + LOG(SUM(posts))

		CASE('c')
			
			scidj = ids(cid,4,:)
			CALL uid(scidj, kmax)
			IF (SUM(scidj) > 0) THEN
			  k = 1
			  DO WHILE(scidj(k) /= 0)
				  IF (verbose) THEN 
					WRITE(1,'(A)') 'tracer : sibchild(case="c") calling posterior()'
				  END IF
				  CALL posterior(cid, scidj(k)) !modify pp, ppid and ppk
				  p = pp(ppk,:)
				  posts = posts * (p/SUM(p))
				  weight = weight + LOG(SUM(p))
				  k = k + 1
			  END DO
			END IF
			likk = lik(cid,:)
			output = (likk/SUM(likk)) * (posts/SUM(posts))
			weight = weight + LOG(SUM(likk)) + LOG(SUM(posts))

		END SELECT

	END SUBROUTINE sibchild

!--- MISC 

	SUBROUTINE debugfile()

	   IMPLICIT NONE
	   INTEGER :: i
	   
	   LOGICAL, PARAMETER :: printped = .TRUE.
	   LOGICAL, PARAMETER :: printlik = .TRUE.   
	   LOGICAL, PARAMETER :: printids = .TRUE.
	   
	   WRITE(1,'(A,3ES11.3)') 'hwpr = ',hwpr
	   
	   !write the pedigree
	   IF(printped) THEN
		 WRITE(1,'(A)') 'ped = '
		 DO i=1,n
				WRITE(1,'(4I4)') ped(i,:)
		 END DO
	   END IF
	   
	   !write the likelihood matrix
	   IF(printlik) THEN
		 WRITE(1,'(A)') 'likelihood = '
		 DO i=1,n
				WRITE(1,'(3ES11.3)') lik(i,:)
		 END DO
	   END IF
	   
	  !write ids table
	  IF(printids) THEN
		 WRITE(1,'(A)') 'ids = '
		 DO i=1,n
				WRITE(1,'(A,i4)')   'indid : ',    ped(i,1) 
				WRITE(1,'(A,i4)')   'motherid : ', ped(i,4) 
				WRITE(1,'(A,i4)')   'fatherid : ', ped(i,3) 
				WRITE(1,*) 'smid : ',   ids(i,1,:)
				WRITE(1,*) 'sfid : ',   ids(i,2,:)
				WRITE(1,*) 'fsibid : ', ids(i,3,:)
				WRITE(1,*) 'sid : ',    ids(i,4,:)
				WRITE(1,*) 'cid : ',    ids(i,5,:)
				WRITE(1,'(A)') '' !saute une ligne
		 END DO
	   END IF
	  
	END SUBROUTINE debugfile
	
	SUBROUTINE vouterprod(a, b, n, m, c)
	 !outer product a(1 ... n) %o% b(1 ... m) vectorized (equiv R >  c <- as.vector(outer(a,b))
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: a
		DOUBLE PRECISION, DIMENSION(m), INTENT(IN) :: b
		DOUBLE PRECISION, DIMENSION(n*m) :: c
		INTEGER :: n, m, k, j
		k = 1
		DO j=1,m
		  c(k:(k+n-1)) = a(:) * b(j)
		  k = k + n
		END DO
	   
	END SUBROUTINE vouterprod
 
	 INTEGER FUNCTION posint(x,table,n)
	   IMPLICIT NONE
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

	LOGICAL FUNCTION inint(x,table,n)
	  IMPLICIT NONE
		! return true if x is in table
		INTEGER, DIMENSION(n), INTENT(IN) :: table
		INTEGER, INTENT(IN) :: x, n
		INTEGER :: pos
		pos = 1
		inint = .FALSE.
		DO WHILE ((.NOT. inint) .AND. (pos <= n))
		  if (table(pos)==x) THEN
			inint = .TRUE.
		  ELSE
			pos = pos + 1
		  END IF
		END DO 
	END FUNCTION inint

	SUBROUTINE rescale(pr, n)
	  IMPLICIT NONE
		INTEGER, DIMENSION(n) :: pr
		INTEGER :: n
		pr = pr / SUM(pr)
	END SUBROUTINE
   
	SUBROUTINE uid(table,n)
	  IMPLICIT NONE
		! return the unique ids of table (remove duplicates)
		! and replace other values by 0 
		INTEGER, DIMENSION(n) :: table,tu
		INTEGER :: i,nu,n
		!LOGICAL :: inint 
		nu = 1
		tu(1) = table(1)
		DO i=2,n
		  IF (.NOT. inint(table(i),tu(1:nu),nu)) THEN
			!found a new unique value 
			nu = nu + 1
			tu(nu) = table(i)
		  END IF
		END DO
		table = tu
		table((nu + 1):n) = 0 
	END SUBROUTINE uid

     
 END SUBROUTINE peeling