
SUBROUTINE genet_f( hwpr, gprv, ngen, af, nl )
   IMPLICIT NONE
   DOUBLE PRECISION hwpr(*)
   DOUBLE PRECISION gprv(*)
   INTEGER ngen(*)
   DOUBLE PRECISION af(*)
   INTEGER nl(*)
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: gpr
   INTEGER :: ng
   
   ng = ngen(1)
   ALLOCATE(gpr(ng,ng,ng))
   CALL compute_HWpr(af(1:ng), nl(1), hwpr(1:ng)) ! compute hwpr
   CALL compute_TRpr(nl(1), gpr) ! compute gpr
   gprv(1:(ng**3)) = PACK(gpr,  mask = .true.)
   DEALLOCATE(gpr)

END SUBROUTINE genet_f
  
  