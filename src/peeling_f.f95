 SUBROUTINE peeling_f( counsid, pedv, likv, counspr, hwpr, gprv, idsv, w, ll, ngen, ni, km)  
  
  IMPLICIT NONE
  INTEGER counsid(*)
  INTEGER pedv(*)
  DOUBLE PRECISION likv(*)
  DOUBLE PRECISION counspr(*)
  DOUBLE PRECISION hwpr(*)
  DOUBLE PRECISION gprv(*)
  INTEGER idsv(*)
  DOUBLE PRECISION w(*)
  DOUBLE PRECISION ll(*)
  INTEGER ngen(*)
  INTEGER ni(*)
  INTEGER km(*)
  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ped
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ids
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: gpr
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: lik
  INTEGER :: ng,n,kmax

  n = ni(1)
  ng = ngen(1)
  kmax = km(1)

  ALLOCATE(ped(n,4))
  ALLOCATE(lik(n,ng))
  ALLOCATE(ids(n,5,kmax))
  ALLOCATE(gpr(ng,ng,ng))

  ped     = RESHAPE(pedv(1:n*4), (/ n,4 /) )
  lik     = RESHAPE(likv(1:n*ng), (/ n,ng /) )
  gpr     = RESHAPE(gprv(1:ng*ng*ng), (/ ng,ng,ng /) )
  ids     = RESHAPE(idsv(1:n*5*kmax), (/ n,5,kmax /) )

  CALL peeling(ll(1), n, ng, kmax, ped, ids, lik, counspr(1:ng), w(1), counsid(1), hwpr(1:ng), gpr) ! modifies counspr, w and ll

  DEALLOCATE(ped)
  DEALLOCATE(lik)
  DEALLOCATE(ids)
  DEALLOCATE(gpr)
  
END SUBROUTINE peeling_f