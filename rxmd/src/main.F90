!------------------------------------------------------------------------------
program rxmd
use atom_vars; use bo; use md_context; use rxmd_params; use cmdline_args; use mpi_vars
use ff_params; use energy_terms; use support_funcs; use comms; use fileio_funcs
use init_funcs; use qeq_vars; use qeq_funcs
!use CG !!FIXME!!
!------------------------------------------------------------------------------
implicit none
integer :: i,it1,it2,irt,provided, nstep, ierr
real(8) :: ctmp, dr(3)

type(md_context_type) :: mcx
type(forcefield_params) :: ffp
type(atom_var_type) :: avs
type(qeq_var_type) :: qvt
type(rxmd_param_type) :: rxp
type(cmdline_arg_type) :: cla
type(mpi_var_type) :: mpt
type(bo_var_type) :: bos

call GetCmdLineArgs(cla)

call GetMPIVariables(mpt)

call initialize_md_context(mcx)

if(mpt%myid==0)  print'(a30)', 'rxmd has started'

!--- read ffield file
CALL GETPARAMS(ffp, cla%FFPath, mcx%FFDescript)

!--- initialize the MD system
CALL INITSYSTEM(mcx, ffp, avs, qvt, bos, rxp, cla, mpt)

!!FIXME!!
!if(rxp%mdmode==10) call ConjugateGradient(ffp, mpt, bos, avs, avs%atype, avs%pos, rxp%ftol)

call QEq(ffp, avs, qvt, mpt, rxp, mcx)
call FORCE(ffp, mpt, bos, avs, qvt, mcx)

!--- Enter Main MD loop 
call system_clock(it1,irt)

do nstep=0, rxp%ntime_step-1

!--- FIXME: user-defined type var cannot be used as a loop counter. update time of md context here. 
   mcx%nstep=nstep

   if(mod(nstep, rxp%pstep)==0) then
       call PRINTE(mcx, mpt, rxp%pstep, avs%atype, avs%v, avs%q)
   endif

   if(mod(nstep, rxp%fstep)==0) &
        call OUTPUT(mcx, ffp, avs, qvt, bos, rxp, mpt, GetFileNameBase(cla%dataDir, mcx%current_step+nstep))

   if(mod(nstep, rxp%sstep)==0 .and. rxp%mdmode==4) &
      avs%v(1:mcx%NATOMS,1:3)=rxp%vsfact*avs%v(1:mcx%NATOMS,1:3)

   if(mod(nstep,rxp%sstep)==0 .and. rxp%mdmode==5) then
      ctmp = (rxp%treq*UTEMP0)/( mcx%GKE*UTEMP )
      avs%v(1:mcx%NATOMS,1:3)=sqrt(ctmp)*avs%v(1:mcx%NATOMS,1:3)
   endif

   if(mod(nstep,rxp%sstep)==0.and.(rxp%mdmode==0.or.rxp%mdmode==6)) &
      call INITVELOCITY(mcx, ffp, avs%atype, avs%v, rxp%treq)

   if(mod(nstep,rxp%sstep)==0.and.rxp%mdmode==7) &
      call ScaleTemperature(mcx, ffp, mpt, rxp%treq, avs%atype, avs%v)

!--- update velocity
   call vkick(mcx, 1.d0, avs%atype, avs%v, avs%f) 

!--- update coordinates
   qvt%qsfv(1:mcx%NATOMS) = qvt%qsfv(1:mcx%NATOMS) + &
      0.5d0 * rxp%dt * qvt%Lex_w2 * (avs%q(1:mcx%NATOMS)-qvt%qsfp(1:mcx%NATOMS))

   qvt%qsfp(1:mcx%NATOMS) = qvt%qsfp(1:mcx%NATOMS) + rxp%dt * qvt%qsfv(1:mcx%NATOMS)

   avs%pos(1:mcx%NATOMS,1:3) = avs%pos(1:mcx%NATOMS,1:3) + rxp%dt * avs%v(1:mcx%NATOMS,1:3)

!--- migrate atoms after positions are updated
   call COPYATOMS(mcx, avs, qvt, mpt, MODE_MOVE, [0.d0, 0.d0, 0.d0])
   
   if(mod(nstep,rxp%qstep)==0) call QEq(ffp, avs, qvt, mpt, rxp, mcx)
   call FORCE(ffp, mpt, bos, avs, qvt, mcx)

!--- update velocity
   call vkick(mcx, 1.d0, avs%atype, avs%v, avs%f) 

   qvt%qsfv(1:mcx%NATOMS) = qvt%qsfv(1:mcx%NATOMS) +  &
      0.5d0 * rxp%dt * qvt%Lex_w2 * (avs%q(1:mcx%NATOMS)-qvt%qsfp(1:mcx%NATOMS))

enddo

!--- FIXME: user-defined type var cannot be used as a loop counter. update time of md context here. 
mcx%nstep=nstep

!--- save the final configurations
call OUTPUT(mcx, ffp, avs, qvt, bos, rxp, mpt, GetFileNameBase(cla%dataDir, mcx%current_step+mcx%nstep))

!--- update rxff.bin in working directory for continuation run
if(rxp%isBinary) call WriteBIN(mcx, avs, qvt, rxp, mpt, GetFileNameBase(cla%dataDir, -1))

call system_clock(it2,irt)
mcx%it_timer(Ntimer)=(it2-it1)

call FinalizeMD(mcx, mpt%myid, irt)

call MPI_FINALIZE(ierr)
end PROGRAM

!------------------------------------------------------------------------------
subroutine FinalizeMD(mcx, myid, irt)
use md_context; use MemoryAllocator; use mpi_vars
!------------------------------------------------------------------------------
implicit none

type(md_context_type),intent(in) :: mcx
integer,intent(in) :: myid
integer,intent(in) :: irt ! time resolution
integer,allocatable :: ibuf(:),ibuf1(:)

integer :: it_timer_max(Ntimer), it_timer_min(Ntimer)

integer :: i, ierr

allocate(ibuf(nmaxas),ibuf1(nmaxas))
ibuf(:)=0
do i=1,nmaxas
   ibuf(i)=maxval(mcx%maxas(:,i))
enddo
call MPI_ALLREDUCE(ibuf, ibuf1, nmaxas, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

call MPI_ALLREDUCE(mcx%it_timer, it_timer_max, Ntimer, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(mcx%it_timer, it_timer_min, Ntimer, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

if(myid==0) then
   print'(a)','----------------------------------------------'
   print'(a20,i12)', 'MAXNEIGHBS: ', ibuf1(2)
   print'(a20,i12)', 'MAXNEIGHBS10: ', ibuf1(3)
   print'(a20,i12)', 'MAXNBUFFER(MOVE): ', ibuf1(1)+ibuf1(4)
   print'(a20,i12)', 'MAXNBUFFER(COPY): ', ibuf1(1)+ibuf1(5)
   print'(a20,i12)', 'QEq Iterations: ', it_timer_max(24)
   print'(a20,f12.2)','Memory (MB): ', GetTotalMemory()*1d-6
   print*

   print'(a20,f12.4,3x,f12.4)','QEq: ',  dble(it_timer_max(1))/irt, dble(it_timer_min(1))/irt
   print'(a20,f12.4,3x,f12.4)','qeq_initialize: ',  dble(it_timer_max(16))/irt, dble(it_timer_min(16))/irt
   print'(a20,f12.4,3x,f12.4)','get_hsh: ',  dble(it_timer_max(18))/irt, dble(it_timer_min(18))/irt
   print'(a20,f12.4,3x,f12.4)','get_gradient: ',  dble(it_timer_max(19))/irt, dble(it_timer_min(19))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','LINKEDLIST: ',  dble(it_timer_max(3))/irt, dble(it_timer_min(3))/irt
   print'(a20,f12.4,3x,f12.4)','COPYATOMS: ',    dble(it_timer_max(4))/irt, dble(it_timer_min(4))/irt
   print'(a20,f12.4,3x,f12.4)','send_rec: ', dble(it_timer_max(25))/irt, dble(it_timer_min(25))/irt
   print'(a20,f12.4,3x,f12.4)','store_atoms: ', dble(it_timer_max(26))/irt, dble(it_timer_min(26))/irt
   print'(a20,f12.4,3x,f12.4)','append_atoms: ', dble(it_timer_max(27))/irt, dble(it_timer_min(27))/irt

   print'(a20,f12.4,3x,f12.4)','NEIGHBORLIST: ', dble(it_timer_max(5))/irt, dble(it_timer_min(5))/irt
   print'(a20,f12.4,3x,f12.4)','GetNBPairList: ', dble(it_timer_max(15))/irt, dble(it_timer_min(15))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','BOCALC: ', dble(it_timer_max(6))/irt, dble(it_timer_min(6))/irt
   print'(a20,f12.4,3x,f12.4)','ENbond: ', dble(it_timer_max(7))/irt, dble(it_timer_min(7))/irt
   print'(a20,f12.4,3x,f12.4)','Ebond: ', dble(it_timer_max(8))/irt, dble(it_timer_min(8))/irt
   print'(a20,f12.4,3x,f12.4)','Elnpr: ', dble(it_timer_max(9))/irt, dble(it_timer_min(9))/irt
   print'(a20,f12.4,3x,f12.4)','Ehb: ', dble(it_timer_max(10))/irt, dble(it_timer_min(10))/irt
   print'(a20,f12.4,3x,f12.4)','E3b: ', dble(it_timer_max(11))/irt, dble(it_timer_min(11))/irt
   print'(a20,f12.4,3x,f12.4)','E4b: ', dble(it_timer_max(12))/irt, dble(it_timer_min(12))/irt
   print'(a20,f12.4,3x,f12.4)','ForceBondedTerms: ', dble(it_timer_max(13))/irt, dble(it_timer_min(13))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','WriteBND: ', dble(it_timer_max(20))/irt, dble(it_timer_min(20))/irt
   print'(a20,f12.4,3x,f12.4)','WritePDB: ', dble(it_timer_max(21))/irt, dble(it_timer_min(21))/irt
   print'(a20,f12.4,3x,f12.4)','ReadBIN: ', dble(it_timer_max(22))/irt, dble(it_timer_min(22))/irt
   print'(a20,f12.4,3x,f12.4)','WriteBIN: ', dble(it_timer_max(23))/irt, dble(it_timer_min(23))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','total (sec): ',dble(it_timer_max(Ntimer))/irt, dble(it_timer_min(Ntimer))/irt

   print'(a)','----------------------------------------------'

   print'(a30)', 'rxmd successfully finished'
endif

deallocate(ibuf,ibuf1)

end subroutine
