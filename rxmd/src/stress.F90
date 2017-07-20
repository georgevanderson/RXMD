!--------------------------------------------------------------------------------------------------------------
subroutine stress(mcx, mass, atype, pos, v, f)
use md_context; use ff_params; use mpi_vars
! calculate stress tensor components of kinetic energy part. 
!--------------------------------------------------------------------------------------------------------------
implicit none

real(8),allocatable,dimension(:) :: mass

type(md_context_type) :: mcx
real(8),allocatable :: atype(:),pos(:,:),v(:,:),f(:,:)

integer :: i,ity, j,jty,icmp, m,n, ierr
integer :: c1,c2,c3,c4,c5,c6
real(8) :: dr(3),dr2
real(8) :: mh, cstr(0:6,10)
real(8) :: dbuf(6), Gdbuf(6)
!--- cutoff length of stress calculation range
real(8) :: rcstr2 = 5.d0**2
real(8) :: qdummy(1)

!--- update buffer atom's stress value
dr=4*mcx%lcsize(1:3)
call COPYATOMS(mcx, atype, pos, v, f, qdummy, MODE_STRESSCALC, dr)

!--- get the potential contribution of internal pressure 
do i=1, mcx%NATOMS
  ity = atype(i)
  mh = mass(ity)
  mcx%xx = mh*v(i,1)*v(i,1)
  mcx%yy = mh*v(i,2)*v(i,2)
  mcx%zz = mh*v(i,3)*v(i,3)
  mcx%yz = mh*v(i,2)*v(i,3)
  mcx%zx = mh*v(i,1)*v(i,3)
  mcx%xy = mh*v(i,1)*v(i,2)

!--- one step values.
  mcx%astr(1,i) = mcx%astr(1,i) + mcx%xx
  mcx%astr(2,i) = mcx%astr(2,i) + mcx%yy
  mcx%astr(3,i) = mcx%astr(3,i) + mcx%zz
  mcx%astr(4,i) = mcx%astr(4,i) + mcx%yz
  mcx%astr(5,i) = mcx%astr(5,i) + mcx%zx
  mcx%astr(6,i) = mcx%astr(6,i) + mcx%xy
enddo

dbuf(1:6)=0.d0
do i=1, mcx%NATOMS
   dbuf(1:6)=dbuf(1:6)+mcx%astr(1:6,i)
enddo

call MPI_ALLREDUCE(dbuf,Gdbuf,size(dbuf),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

dbuf(1:6)=Gdbuf(1:6)/mcx%MDBOX
mcx%pint(1,1)=dbuf(1) !xx
mcx%pint(2,2)=dbuf(2) !yy
mcx%pint(3,3)=dbuf(3) !zz
mcx%pint(3,1)=dbuf(4) !zx
mcx%pint(1,3)=dbuf(4) !xz
mcx%pint(2,3)=dbuf(5) !yz
mcx%pint(3,2)=dbuf(5) !zy
mcx%pint(1,2)=dbuf(6) !xy
mcx%pint(2,1)=dbuf(6) !yx

return
end subroutine
!--------------------------------------------------------------------------------------------------------------
