
SUBROUTINE preproc_f(pedbrutv , pedv, idsv, disv, agev, tesv, ascv, ngen, ni, ndi, nl, km)

  IMPLICIT NONE
  
  INTEGER pedbrutv(*)
  INTEGER pedv(*)
  INTEGER idsv(*)
  INTEGER disv(*)
  INTEGER agev(*)
  DOUBLE PRECISION tesv(*)
  INTEGER ascv(*)
  INTEGER ngen(*)
  INTEGER ni(*)
  INTEGER ndi(*)
  INTEGER nl(*)
  INTEGER km(*)
  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: pedbrut, ped, tes, dis, asc, age
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ids
  INTEGER :: ng,n,nloci,ndis,kmax

  n = ni(1)
  ng = ngen(1)
  kmax = km(1)
  ndis = ndi(1)
  nloci = nl(1)

  ALLOCATE(ped(n,4))
  ALLOCATE(pedbrut(n,5+nloci+2*ndis))
  ALLOCATE(tes(n,ng))
  ALLOCATE(dis(n,ndis))
  ALLOCATE(age(n,ndis))
  ALLOCATE(asc(n,ng))
  ALLOCATE(ids(n,5,kmax))

  pedbrut = RESHAPE(pedbrutv(1:n*(5+nloci+2*ndis)), (/ n, 5+nloci+2*ndis /) )

  CALL preproc(n,ng,nloci,ndis,pedbrut,ped,kmax,ids,dis,tes,asc,age) !modifies ped, ids, dis, tes, asc

  idsv(1:n*5*kmax) = PACK(ids, mask = .true.)
  pedv(1:4*n)      = PACK(ped, mask = .true.)
  agev(1:n*ndis)   = PACK(age, mask = .true.)
  tesv(1:n*ng)     = PACK(tes, mask = .true.)
  ascv(1:n*ng)     = PACK(asc, mask = .true.)
  disv(1:n*ndis)   = PACK(dis, mask = .true.)

  DEALLOCATE(ped)
  DEALLOCATE(pedbrut)
  DEALLOCATE(ids)
  DEALLOCATE(dis)
  DEALLOCATE(age)
  DEALLOCATE(tes)
  DEALLOCATE(asc)

END SUBROUTINE preproc_f
