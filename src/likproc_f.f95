SUBROUTINE likproc_f( pedbrutv, ftv, disv, tesv, ascv, agev, liknumv, likdenomv, ngen, ni, ndi, nl )
  
  IMPLICIT NONE
  
  INTEGER pedbrutv(*)
  DOUBLE PRECISION ftv(*)
  INTEGER disv(*)
  INTEGER tesv(*)
  INTEGER ascv(*)
  INTEGER agev(*)
  DOUBLE PRECISION liknumv(*)
  DOUBLE PRECISION likdenomv(*)
  INTEGER ngen(*)
  INTEGER ni(*)
  INTEGER ndi(*)
  INTEGER nl(*)

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: pedbrut, tes, dis, asc, age
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: ft
  INTEGER :: ng, n, nloci, ndis
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: liknum, likdenom

  n = ni(1)
  ng = ngen(1)
  ndis = ndi(1)
  nloci = nl(1)

  ALLOCATE(pedbrut(n,5+nloci+2*ndis))
  ALLOCATE(tes(n,ng))
  ALLOCATE(dis(n,ndis))
  ALLOCATE(age(n,ndis))
  ALLOCATE(asc(n,ng))
  ALLOCATE(ft(121,2,ndis,ng))
  ALLOCATE(liknum(n,ng))
  ALLOCATE(likdenom(n,ng))

  pedbrut = RESHAPE(pedbrutv(1:n*(5+nloci+2*ndis)), (/ n, 5+nloci+2*ndis /) )
  ft      = RESHAPE(ftv(1:121*2*ndis*ng), (/ 121,2,ndis,ng /) )
  age     = RESHAPE(agev(1:n*ndis), (/ n,ndis /) )
  dis     = RESHAPE(disv(1:n*ndis), (/ n,ndis /) )
  tes     = RESHAPE(tesv(1:n*ng), (/ n,ng /) )
  asc     = RESHAPE(ascv(1:n*ng), (/ n,ng /) )

  CALL likproc(pedbrut, n, ng, ndis, nloci, ft, dis, age, tes, asc, liknum, likdenom)

  liknumv(1:n*ng)   = PACK(liknum, mask = .true.)
  likdenomv(1:n*ng) = PACK(likdenom, mask = .true.)

  DEALLOCATE(liknum)
  DEALLOCATE(likdenom)
  DEALLOCATE(pedbrut)
  DEALLOCATE(dis)
  DEALLOCATE(age)
  DEALLOCATE(tes)
  DEALLOCATE(asc)
  DEALLOCATE(ft)
  
END SUBROUTINE likproc_f
