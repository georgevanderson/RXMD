module list_funcs

contains
!----------------------------------------------------------------------------------------
subroutine LINKEDLIST(mcx, atype, rreal, cellDims, headAtom, atomList, NatomPerCell, Ncells, NLAYERS)
use md_context; use support_funcs
! partitions the volume into linked-list cells <lcsize>
!----------------------------------------------------------------------------------------
implicit none

type(md_context_type),intent(inout) :: mcx
real(8),allocatable,intent(in) :: atype(:), rreal(:,:)
real(8),intent(in) :: cellDims(3)

integer,intent(in) :: Ncells(3), NLAYERS
integer,intent(inout),allocatable :: atomList(:)
integer,intent(inout) :: NatomPerCell(-NLAYERS:Ncells(1)-1+NLAYERS, &
                                    -NLAYERS:Ncells(2)-1+NLAYERS, &
                                    -NLAYERS:Ncells(3)-1+NLAYERS) 
integer,intent(inout) :: headAtom(-NLAYERS:Ncells(1)-1+NLAYERS, & 
                                -NLAYERS:Ncells(2)-1+NLAYERS, &
                                -NLAYERS:Ncells(3)-1+NLAYERS) 

real(8),allocatable :: rnorm(:,:)
integer :: n, l(3), j

integer :: ti,tj,tk
call system_clock(ti,tk)

if(.not.allocated(rnorm) ) allocate(rnorm(mcx%NBUFFER,3))

call xu2xs(mcx, rreal,rnorm,mcx%copyptr(6))

headAtom(:,:,:) = -1; atomList(:) = 0; NatomPerCell(:,:,:)=0

!--- copyptr(6) stores the last atom index copied in COPYATOMS.
do n=1, mcx%copyptr(6) 

   if(nint(atype(n))==0) cycle

   l(1:3) = floor(rnorm(n,1:3)/cellDims(1:3))

   atomList(n) = headAtom(l(1), l(2), l(3))
   headAtom(l(1), l(2), l(3)) = n
   NatomPerCell(l(1), l(2), l(3)) = NatomPerCell(l(1), l(2), l(3)) + 1
enddo

call system_clock(tj,tk)
mcx%it_timer(3)=mcx%it_timer(3)+(tj-ti)

end subroutine 

!----------------------------------------------------------------------
subroutine NEIGHBORLIST(mcx, ffp, nlayer, atype, pos)
use md_context; use ff_params
! calculate neighbor list for atoms witin cc(1:3, -nlayer:nlayer) cells.
!----------------------------------------------------------------------
implicit none

type(md_context_type),intent(inout) :: mcx
type(forcefield_params),intent(in) :: ffp
integer,intent(in) :: nlayer
real(8),intent(in),allocatable :: atype(:), pos(:,:)

integer :: c1,c2,c3, ic(3), c4, c5, c6
integer :: n, n1, m, m1, nty, mty, inxn
real(8) :: dr(3), dr2

integer :: i,j,i1,j1, ierr
logical :: isFound

integer :: ti,tj,tk
call system_clock(ti,tk)

mcx%nbrlist(:,0) = 0

!$omp parallel do default(shared) collapse(3) & 
!$omp private(c1,c2,c3,ic,c4,c5,c6,n,n1,m,m1,nty,mty,inxn,dr,dr2) 
DO c1=-nlayer, mcx%cc(1)-1+nlayer
DO c2=-nlayer, mcx%cc(2)-1+nlayer
DO c3=-nlayer, mcx%cc(3)-1+nlayer

  m = mcx%header(c1, c2, c3)
  do m1=1, mcx%nacell(c1, c2, c3)
     mty = nint(atype(m))

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = [c1, c2, c3] + [c4, c5, c6]

        n = mcx%header(ic(1),ic(2),ic(3))
        do n1=1, mcx%nacell(ic(1), ic(2), ic(3))

           if(n/=m) then
             nty = nint(atype(n))
             inxn = ffp%inxn2(mty, nty)

             dr(1:3) = pos(n,1:3) - pos(m,1:3) 
             dr2 = sum(dr(1:3)*dr(1:3))

             if(dr2<mcx%rc2(inxn)) then 
                mcx%nbrlist(m, 0) = mcx%nbrlist(m, 0) + 1
                mcx%nbrlist(m, mcx%nbrlist(m, 0)) = n
             endif 
           endif

           n=mcx%llist(n) 
        enddo
     enddo; enddo; enddo

     m=mcx%llist(m)
  enddo
enddo; enddo; enddo
!$omp end parallel do 

!--- to get the reverse information (i.e. from i,j1&j to i1), store <i1> into <nbrindx>.

!$omp parallel do default(shared) private(i,i1,j,j1,isFound)
do i=1, mcx%copyptr(6)
   do i1 = 1, mcx%nbrlist(i,0)
      j = mcx%nbrlist(i,i1)
      isFound=.false.
      do j1 = 1, mcx%nbrlist(j,0)
         if(i == mcx%nbrlist(j,j1)) then
            mcx%nbrindx(i,i1)=j1
            isFound=.true.
         endif
      enddo
      if(.not.isFound) &

      !print'(a,i6,30i4)','ERROR: inconsistency between nbrlist and nbrindx found', &
      !     myid, i,mcx%nbrlist(i,0:mcx%nbrlist(i,0)), j, mcx%nbrlist(j,0:mcx%nbrlist(j,0))
      !! FIXME !! this function is called from FORCE() only and all functions under FORCE() doesn't need to use MPI. 
      ! how do we print myid without passing mpi_var_type? 
      print'(a,30i4)','ERROR: inconsistency between nbrlist and nbrindx found', &
           i,mcx%nbrlist(i,0:mcx%nbrlist(i,0)), j, mcx%nbrlist(j,0:mcx%nbrlist(j,0))
   enddo
enddo
!$omp end parallel do

!--- error trap
n=maxval(mcx%nbrlist(1:mcx%NATOMS,0))
if(n > MAXNEIGHBS) then
   !! FIXME !! this function is called from FORCE() only and all functions under FORCE() doesn't need to use MPI. 
   ! how do we print myid without passing mpi_var_type? 
   !write(6,'(a45,2i5)') "ERROR: overflow of max # in neighbor list, ", myid, n
   write(6,'(a45,i5)') "ERROR: overflow of max # in neighbor list, ", n
   call MPI_FINALIZE(ierr)
   stop
endif

!! FIXME !! this function is called from FORCE() only and all functions under FORCE() doesn't need to use MPI. 
! how do we print myid without passing mpi_var_type? 
!!--- for array size stat
!if(mod(nstep,rxp%pstep)==0) then
!  maxas(nstep/pstep+1,2)=maxval(mcx%nbrlist(1:mcx%NATOMS,0))
!endif

call system_clock(tj,tk)
mcx%it_timer(5)=mcx%it_timer(5)+(tj-ti)

end subroutine

!----------------------------------------------------------------------
subroutine GetNonbondingPairList(mcx, rctap2, pos)
use md_context; use ff_params 
!----------------------------------------------------------------------
implicit none

type(md_context_type),intent(inout) :: mcx
real(8),intent(in) :: rctap2
real(8),intent(in),allocatable :: pos(:,:)

integer :: c1,c2,c3,c4,c5,c6,i,j,m,n,mn,iid,jid
integer :: l2g
real(8) :: dr(3), dr2

integer :: ti,tj,tk
call system_clock(ti,tk)

! reset non-bonding pair list
mcx%nbplist(:,0)=0

!$omp parallel do default(shared),private(c1,c2,c3,c4,c5,c6,i,j,m,n,mn,iid,jid,dr,dr2)
do c1=0, mcx%nbcc(1)-1
do c2=0, mcx%nbcc(2)-1
do c3=0, mcx%nbcc(3)-1

   i = mcx%nbheader(c1,c2,c3)
   do m = 1, mcx%nbnacell(c1,c2,c3)

      do mn = 1, mcx%nbnmesh
         c4 = c1 + mcx%nbmesh(1,mn)
         c5 = c2 + mcx%nbmesh(2,mn)
         c6 = c3 + mcx%nbmesh(3,mn)

         j = mcx%nbheader(c4,c5,c6)
         do n=1, mcx%nbnacell(c4,c5,c6)

            !if(i<j .or. NATOMS<j) then
            if(i/=j) then
               dr(1:3) = pos(i,1:3) - pos(j,1:3)
               dr2 = sum(dr(1:3)*dr(1:3))

               if(dr2<=rctap2) then
                 mcx%nbplist(i,0)=mcx%nbplist(i,0)+1
                 mcx%nbplist(i,mcx%nbplist(i,0))=j
               endif

            endif

            j=mcx%nbllist(j)
         enddo
       enddo

      i=mcx%nbllist(i)
   enddo
enddo; enddo; enddo
!$omp end parallel do

call system_clock(tj,tk)
mcx%it_timer(15)=mcx%it_timer(15)+(tj-ti)

end subroutine

end module 
