!---------------------------------------------------------------------------------
module CG
use atom_vars; use energy_terms; use bo; use qeq_vars; use qeq_funcs; use comms
!---------------------------------------------------------------------------------
integer,parameter :: CG_MaxMinLoop = 500
integer,parameter :: CG_MaxLineMinLoop = 100
integer,parameter :: CG_MaxBracketLoop = 20
integer,parameter :: CG_MaxGSLoop = 100


!--- Wolfe condition parameters
real(8),parameter :: CG_WC1 = 1d-4
real(8),parameter :: CG_WC2 = 0.1d0

!--- golden section tolerance
real(8),parameter :: CG_GStol = 1d-6

!--- conjugate gradient tolerance. ftol in rxmd.in
real(8) :: CG_tol = 1d-4

real(8),parameter :: CG_EPS= 1d-16 ! a check to emit warning message

!--- these vars are mostly used as placeholder and allocated in the CG initializer.
type(atom_var_type) :: avsTmp
type(bo_var_type) :: bos_cg
type(qeq_var_type) :: qvt_cg

contains

!---------------------------------------------------------------------------------
subroutine initialize_cg(mcx)
use md_context; use atom_vars; use bo; use qeq_vars; use qeq_funcs
!---------------------------------------------------------------------------------
   implicit none
   type(md_context_type),intent(inout) :: mcx

   call initialize_atom_vars(avsTmp,mcx%NBUFFER) 
   call initialize_bo_vars(bos_cg, mcx%NBUFFER)
   call initialize_qeq_vars(qvt_cg, mcx%NBUFFER, MAXNEIGHBS10)
   
end subroutine

!---------------------------------------------------------------------------------
subroutine ConjugateGradient(ffp, avs, mcx, rxp, atype, pos, cla, mpt, ftol)
use md_context; use ff_params; use bo; use cmdline_args; use mpi_vars; 
use qeq_vars; use qeq_funcs; use rxmd_params; use ff_params; use fileio_funcs
!---------------------------------------------------------------------------------
implicit none

type(cmdline_arg_type),intent(in) :: cla
type(mpi_var_type),intent(in) :: mpt

type(forcefield_params),intent(in) :: ffp
type(atom_var_type),intent(inout) :: avs
type(md_context_type),intent(inout) :: mcx
type(rxmd_param_type),intent(in) :: rxp

real(8) :: ftol

!real(8) :: atype(NBUFFER),pos(NBUFFER,3)
!real(8) :: f(NBUFFER,3),v(NBUFFER,3),q(NBUFFER)
real(8),allocatable :: atype(:),pos(:,:)
real(8),allocatable :: f(:,:),v(:,:),q(:)

real(8) :: p(mcx%NBUFFER,3) ! search direction
real(8) :: gold(mcx%NBUFFER,3),gnew(mcx%NBUFFER,3) ! old and new gradients
real(8) :: GPE(0:13),GPEold,GPEnew

real(8) :: vnorm, stepl, beta1,beta2,beta3

integer :: cgLoop, i

CG_tol = ftol
v(:,:)=0.d0

if(mpt%myid==0) print'(a40,1x,es10.2)', NEW_LINE('A')//'Start structural optimization.', ftol


call QEq(ffp, avs, qvt_cg, mpt, rxp, mcx)
call FORCE(ffp, mpt, bos_cg, avs, qvt_cg, mcx)
gnew(1:mcx%NATOMS,1:3) = avs%f(1:mcx%NATOMS,1:3) !!!!! FIXME, need to return gnew from FORCE()

!--- initialize search direction with gradient
p(1:mcx%NATOMS,1:3)=gnew(1:mcx%NATOMS,1:3)

mcx%PE(0)=sum(mcx%PE(1:13))
call MPI_ALLREDUCE(mcx%PE, GPE, size(mcx%PE), MPI_DOUBLE_PRECISION, MPI_SUM, mpt%mycomm, mpt%ierr)
GPEnew=GPE(0)

!--- if no bracket range was found here, you are at the energy minimum already. 
!call BracketSearchRange(mpt%myid,avs%atype,avs%pos,p,stepl)
call BracketSearchRange(avs,ffp,mcx,rxp,mpt, avs%atype,avs%pos,p,stepl)

do cgLoop = 0, CG_MaxMinLoop-1

   !call LineMinimization(mpt%myid,avs%atype,avs%pos,p,gnew,stepl)
   call LineMinimization(avs,ffp,mcx,rxp,mpt, avs%atype,avs%pos,p,gnew,stepl)
   gold(1:mcx%NATOMS,1:3)=gnew(1:mcx%NATOMS,1:3)

   call QEq(ffp, avs, qvt_cg, mpt, rxp, mcx)
   call FORCE(ffp, mpt, bos_cg, avs, qvt_cg, mcx)

   call OUTPUT(mcx, ffp, avs, qvt_cg, bos_cg, rxp, mpt, GetFileNameBase(cla%dataDir, cgLoop))

   GPEold=GPEnew
   mcx%PE(0)=sum(mcx%PE(1:13))
   call MPI_ALLREDUCE(mcx%PE, GPE, size(mcx%PE), MPI_DOUBLE_PRECISION, MPI_SUM, mpt%mycomm, mpt%ierr)
   GPEnew=GPE(0)

   if(abs(GPEnew-GPEold)<=CG_tol*mcx%GNATOMS) then
      if(mpt%myid==0) print'(a30,i6)','Energy minimum was found. Exiting CG loop.', cgLoop

      !call OUTPUT(atype, pos, v, q, GetFileNameBase(cla%dataDir, cgLoop))
      call OUTPUT(mcx, ffp, avs, qvt_cg, bos_cg, rxp, mpt, GetFileNameBase(cla%dataDir, cgLoop))
      exit
   endif

   beta1=DotProductVec3D(gold,gold,mcx%NATOMS, mpt)
   beta2=DotProductVec3D(gnew,gnew,mcx%NATOMS, mpt)
   beta3=DotProductVec3D(gnew,gold,mcx%NATOMS, mpt)

   if(mpt%myid==0) print'(a30,i6,4es15.5,2es20.10)', &
     'b1,b2,b3,(b2-b3)/b1: ',cgLoop,beta1,beta2,beta3,(beta2-beta3)/beta1,GPEnew,GPEold

   p(1:mcx%NATOMS,1:3) = (beta2-beta3)/beta1*p(1:mcx%NATOMS,1:3) + gnew(1:mcx%NATOMS,1:3)

   !call BracketSearchRange(ffp,mcx,rxp,mpt,atype,pos,p,stepl)
   call BracketSearchRange(avs,ffp,mcx,rxp,mpt, avs%atype,avs%pos,p,stepl)

enddo 

call MPI_Finalize(mpt%ierr)
stop 'successfully finished structural optimization. '

end subroutine ConjugateGradient

!---------------------------------------------------------------------------------
subroutine BracketSearchRange(avs,ffp,mcx,rxp,mpt ,atype,pos,p,stepl) 
use atom_vars; use md_context; use ff_params;use mpi_vars; use rxmd_params;
use fileio_funcs;
! input: atom type, initial coordinate, and search direction.
! output: step length to bracket an energy minimum along the search direction.
!---------------------------------------------------------------------------------
implicit none
!integer,intent(in) :: myid
type(atom_var_type),intent(inout) :: avs
type(forcefield_params),intent(in) :: ffp
type(md_context_type),intent(inout) :: mcx
type(rxmd_param_type),intent(in) :: rxp
type(mpi_var_type),intent(in) :: mpt
!real(8),intent(in) :: atype(NBUFFER),pos(NBUFFER,3),p(NBUFFER,3)
real(8),allocatable,intent(in) :: atype(:),pos(:,:)
real(8),intent(in) :: p(mcx%NBUFFER,3)
real(8),intent(inout) :: stepl
real(8) :: vdummy(mcx%NBUFFER,3), qdummy(mcx%NBUFFER)
integer :: bracketingLoop
logical :: Elower, WolfeC1, WolfeC2

real(8) :: PE

if(mpt%myid==0) print'(a40)', NEW_LINE('A')//'Start BracketSearchRange()'

stepl=1d-2/mcx%GNATOMS; Elower=.true.; WolfeC1=.true.; WolfeC2=.true.

do bracketingLoop = 0, CG_MaxBracketLoop-1

   stepl = stepl*2

   PE = EvaluateEnergyWithStep(ffp,mcx,rxp,mpt, atype,pos,p,stepl)
   call WolfeConditions(avs,ffp,mcx,mpt,rxp, atype,pos,p,stepl,Elower,WolfeC1,WolfeC2)

   if(mpt%myid==0) print'(a30,es15.5,es25.15, 3l3)', &
      'stepl,PE,Elow,Wolfe1,Wolfe2: ', stepl, PE, Elower, WolfeC1, WolfeC2

   if(.not.WolfeC1 .or. .not.WolfeC1) then
      if(mpt%myid==0) print'(a30,es15.5,a1)', 'bracket has been found: ', stepl
      return 
   endif

enddo 

if(mpt%myid==0) print'(a)', 'bracket was not found. terminating the structural optimization'
!call OUTPUT(atype, pos, vdummy, qdummy, "nobracket")

stop

end subroutine BracketSearchRange

!---------------------------------------------------------------------------------
subroutine WolfeConditions(avs,ffp,mcx,mpt,rxp, atype,pos,p,stepl,isLowerEnergy,isArmijoRule,isCurvature)
use atom_vars; use mpi_vars; use ff_params; use md_context; use mpi_vars; use rxmd_params;
! input: atom type, position, search direction and step length
! output: two bools for the Wolfe conditions
! TODO: This function is very similar to EvaluateEnergyWithStep because the function 
! doesn't return force and new NATOM. Can be simplified. 
!---------------------------------------------------------------------------------
implicit none

type(atom_var_type),intent(inout) :: avs
type(forcefield_params),intent(in) :: ffp
type(md_context_type),intent(inout) :: mcx
type(mpi_var_type),intent(in) :: mpt
type(rxmd_param_type),intent(in) :: rxp

real(8),intent(in) :: stepl
real(8),intent(in),allocatable :: atype(:),pos(:,:)
real(8),intent(in) :: p(mcx%NBUFFER,3)
logical,intent(inout) :: isLowerEnergy,isArmijoRule,isCurvature

!real(8) :: atypeTmp(NBUFFER),posTmp(NBUFFER,3),pTmp(NBUFFER,3),qTmp(NBUFFER)
!! FIXME: here we don't really need v but COPYATOM(MOVE) requires it for the move mode. 
!!        thus, I am allocating a 3xNBUFFER dummy array here. a better implementation needed. 
!real(8) :: vdummy(NBUFFER,3) 
real(8) :: GPE(0:13), GPEbefore, GPEafter, fbefore(mcx%NBUFFER,3), fafter(mcx%NBUFFER,3)
real(8) :: pDotdF, pDotdFShift
integer :: NATOMSTmp

! Evaluate df(x) and f(x)
!call QEq(atype, pos, qTmp)
!call FORCE(atype, pos, fbefore, qTmp)
call QEq(ffp, avs, qvt_cg, mpt, rxp, mcx)
call FORCE(ffp, mpt, bos_cg, avs, qvt_cg, mcx)
fbefore(1:mcx%NATOMS,1:3)=avs%f(1:mcx%NATOMS,1:3) !!! FIXME need to return fbefore from FORCE()

mcx%PE(0)=sum(mcx%PE(1:13))
call MPI_ALLREDUCE(mcx%PE, GPE, size(mcx%PE), MPI_DOUBLE_PRECISION, MPI_SUM, mpt%mycomm, mpt%ierr)
GPEbefore=GPE(0)

! Evaluate df(x+alpha*p) and f(x+alpha*p)
avsTmp%pos(1:mcx%NATOMS,1:3)=avs%pos(1:mcx%NATOMS,1:3)+stepl*p(1:mcx%NATOMS,1:3)
avsTmp%atype(1:mcx%NATOMS)=avs%atype(1:mcx%NATOMS)

! FIXME local number of atoms will change after COPYATOM(MOVE). 
! Need to retrieve the origianl value consistent with the coordinates before the move.
NATOMSTmp = mcx%NATOMS

!call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atypeTmp, posTmp, vdummy, fafter, qTmp)
!call QEq(atypeTmp, posTmp, qTmp)
!call FORCE(atypeTmp, posTmp, fafter, qTmp)
call COPYATOMS(mcx, avsTmp, qvt_cg, mpt, MODE_MOVE, [0.d0, 0.d0, 0.d0])
call QEq(ffp, avsTmp, qvt_cg, mpt, rxp, mcx)
call FORCE(ffp, mpt, bos_cg, avsTmp, qvt_cg, mcx)
fafter(1:mcx%NATOMS,1:3)=avsTmp%f(1:mcx%NATOMS,1:3)

! we have p, fbefore, and fafter vectors now. the dot_product between p&fbefore is trivial but 
! p&after requires migration of one of the vectors because there may be atoms that migrate.

mcx%PE(0)=sum(mcx%PE(1:13))
call MPI_ALLREDUCE(mcx%PE, GPE, size(mcx%PE), MPI_DOUBLE_PRECISION, MPI_SUM, mpt%mycomm, mpt%ierr)
GPEafter=GPE(0)

isLowerEnergy = GPEafter < GPEbefore

pDotdF = DotProductVec3D(p,fbefore,mcx%NATOMS, mpt)
isArmijoRule = GPEafter <= GPEbefore + pDotdF*CG_WC1*stepl

! the local number of atoms of f'(x+a) and p may differ after the shift, 
! cannot take their vector dot-product without moving the search vector. 
! piggyback avsTmp%v to migrate p vector here.

mcx%NATOMS = NATOMSTmp 
avsTmp%pos(1:mcx%NATOMS,1:3)=avs%pos(1:mcx%NATOMS,1:3)+stepl*p(1:mcx%NATOMS,1:3)
avsTmp%v(1:mcx%NATOMS,1:3)=p(1:mcx%NATOMS,1:3)
avsTmp%atype(1:mcx%NATOMS)=avs%atype(1:mcx%NATOMS)
!call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0],atypeTmp,posTmp,pTmp,vdummy,qTmp)
call COPYATOMS(mcx, avsTmp, qvt_cg, mpt, MODE_MOVE, [0.d0, 0.d0, 0.d0])

!pDotdFShift = DotProductVec3D(pTmp,fafter,NATOMS)
pDotdFShift = DotProductVec3D(avsTmp%v,fafter,mcx%NATOMS, mpt)
isCurvature = pDotdFShift >= CG_WC2*pDotdF

! recover the original NATOM
mcx%NATOMS = NATOMSTmp 

end subroutine WolfeConditions

!---------------------------------------------------------------------------------
subroutine LineMinimization(avs,ffp,mcx,rxp,mpt, atype,pos,p,g,stepl)
use md_context; use ff_params; use mpi_vars; use rxmd_params;
! Here, we perform line minimization based on golden section search. After obtaining,
! step length, update atom positions and migrate them if they move out of MD box. 
! We also migrate gradient and search vectors according to the position migration 
! in order to perform dot product of old and new gradient vectors. 
! input: atom type, initial coordinate, and search direction.
! output: updated coordinate and associated gradient vector. 
!---------------------------------------------------------------------------------
implicit none

type(atom_var_type),intent(inout) :: avs
type(forcefield_params),intent(in) :: ffp
type(md_context_type),intent(inout) :: mcx
type(rxmd_param_type),intent(in) :: rxp
type(mpi_var_type),intent(in) :: mpt

real(8) :: atype(mcx%NBUFFER),pos(mcx%NBUFFER,3),g(mcx%NBUFFER,3),p(mcx%NBUFFER,3),q(mcx%NBUFFER)
!real(8) :: atypeTmp(NBUFFER),posTmp(NBUFFER,3),fdummy(1,1)
integer :: lineMinLoop, NATOMSTmp
real(8) :: stepl, step0

if(mpt%myid==0) print'(a40)', NEW_LINE('A')//'Start LineMinimization()'

step0=0.d0

call GoldenSectionSearch(ffp,mcx,rxp,mpt, avs%atype,avs%pos,p,step0,stepl)

! First, migrate gradient vector g.
!NATOMStmp = MigrateVec3D(avs%pos,p,g,stepl)
NATOMStmp = MigrateVec3D(mcx, mpt, avs%pos, p, g, stepl)

! Then, migrate atom type, position, and search vector; atype, pos, p.
avs%pos(1:mcx%NATOMS,1:3)=avs%pos(1:mcx%NATOMS,1:3)+stepl*p(1:mcx%NATOMS,1:3)
!call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atype, pos, p, fdummy, q)
call COPYATOMS(mcx, avs, qvt_cg, mpt, MODE_MOVE, [0.d0, 0.d0, 0.d0])

return
end subroutine LineMinimization

!---------------------------------------------------------------------------------
subroutine GoldenSectionSearch(ffp,mcx,rxp,mpt, atype,pos,p,ax,dx)
! ax,dx: left and right boundaries of search range. 
use md_context; use mpi_vars; use ff_params; use rxmd_params;
!---------------------------------------------------------------------------------
implicit none
type(forcefield_params),intent(in) :: ffp
type(md_context_type),intent(inout) :: mcx
type(rxmd_param_type),intent(in) :: rxp
type(mpi_var_type),intent(in) :: mpt

real(8),allocatable,intent(in) :: atype(:),pos(:,:)
real(8),intent(in) :: p(mcx%NBUFFER,3)
real(8) :: ax,bx,cx,dx,PEbx,PEcx
real(8) :: ratio = 1.d0/1.61803398875d0 ! inverse of golden ratio

integer :: GSLoop

if(mpt%myid==0) print'(a30)', 'start golden section step.'

bx=dx-(dx-ax)*ratio
cx=ax+(dx-ax)*ratio
PEbx=EvaluateEnergyWithStep(ffp,mcx,rxp,mpt, atype,pos,p,bx)
PEcx=EvaluateEnergyWithStep(ffp,mcx,rxp,mpt, atype,pos,p,cx)

do GSLoop = 0, CG_MaxLineMinLoop-1

   if(mpt%myid==0) print'(a30,4es15.5,1x,2es25.15)', &
      'ax,bx,cx,dx,PEbx,PEcx: ', ax,bx,cx,dx,PEbx,PEcx

   if(abs(ax-dx)<=CG_GStol/mcx%GNATOMS) then
      if(mpt%myid==0) print'(a30)', 'golden section step finished.'
      exit
   endif

   if(PEbx<PEcx) then
      dx=cx
   else
      ax=bx
   endif

   bx=dx-(dx-ax)*ratio
   cx=ax+(dx-ax)*ratio
   PEbx=EvaluateEnergyWithStep(ffp,mcx,rxp,mpt, atype,pos,p,bx)
   PEcx=EvaluateEnergyWithStep(ffp,mcx,rxp,mpt, atype,pos,p,cx)
enddo

end subroutine GoldenSectionSearch

!---------------------------------------------------------------------------------
subroutine PolynomialFitSearch()
!To be implemented.
!---------------------------------------------------------------------------------
implicit none

end subroutine PolynomialFitSearch

!---------------------------------------------------------------------------------
function MigrateVec3D(mcx, mpt, pos, vec, dir, stepl) result(newNATOMS)
use md_context; use atom_vars; use mpi_vars
!TODO: come up a better way to migrate vectors
!---------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: stepl
type(md_context_type),intent(inout) :: mcx
type(mpi_var_type),intent(in) :: mpt

real(8),allocatable,intent(in) :: pos(:,:)
real(8) :: vec(mcx%NBUFFER,3), dir(mcx%NBUFFER,3)
!real(8) :: atypedummy(mcx%NBUFFER),posTmp(mcx%NBUFFER,3),fdummy(1,1),qdummy(mcx%NBUFFER)
integer :: NATOMSTmp, newNATOMS

!--- keep current NATOMS
NATOMSTmp=mcx%NATOMS

avsTmp%pos(1:mcx%NATOMS,1:3)=pos(1:mcx%NATOMS,1:3)+stepl*dir(1:mcx%NATOMS,1:3)
avsTmp%v(1:mcx%NATOMS,1:3)=vec(1:mcx%NATOMS,1:3)
!call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atypedummy, posTmp, vec, fdummy, qdummy)
call COPYATOMS(mcx, avsTmp, qvt_cg, mpt, MODE_MOVE, [0.d0, 0.d0, 0.d0])

!-- this NATOMS is consistent with the migrated vector.
newNATOMS=mcx%NATOMS

!--- save the original NATOMS that is consistent with input pos. 
mcx%NATOMS=NATOMSTmp

return
end function MigrateVec3D


!---------------------------------------------------------------------------------
function DotProductVec3D(v1, v2, Nelems, mpt) result (ret)
use mpi_vars 
!---------------------------------------------------------------------------------
implicit none

integer,intent(in) :: Nelems
real(8),intent(in) :: v1(Nelems,3), v2(Nelems,3) 
type(mpi_var_type),intent(in) :: mpt

integer :: i
real(8) :: vsum, ret

vsum = 0.d0
do i=1, Nelems
   vsum = vsum + sum(v1(i,1:3)*v2(i,1:3))
enddo

call MPI_ALLREDUCE(vsum, ret, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpt%mycomm, mpt%ierr)

return
end function DotProductVec3D


!---------------------------------------------------------------------------------
subroutine NormalizeVec3D(v, vnorm, nElems, mpt) 
use mpi_vars
!---------------------------------------------------------------------------------
implicit none
type(mpi_var_type),intent(in) :: mpt
real(8),allocatable :: v(:,:)
real(8) :: vsum, vnorm 
integer :: i, nElems

vnorm=sqrt(DotProductVec3D(v,v,nElems,mpt))

if(abs(vnorm)<CG_EPS) then
   if(mpt%myid==0) print'(a,es25.15)', &
    'WARNING: Norm of vector was found too small in NormalizeVector(): vnorm = ', vnorm
endif

v(1:nElems,1:3)=v(1:nElems,1:3)/vnorm

return
end subroutine NormalizeVec3D

!---------------------------------------------------------------------------------
function EvaluateEnergyWithStep(ffp, mcx, rxp, mpt, atype,pos,p,stepl) result(potentialEnergy)
use md_context; use mpi_vars; use ff_params; use rxmd_params;
!---------------------------------------------------------------------------------
implicit none
type(forcefield_params),intent(in) :: ffp
type(rxmd_param_type),intent(in) :: rxp
type(md_context_type),intent(inout) :: mcx
type(mpi_var_type),intent(in) :: mpt

!real(8),intent(in) :: atype(NBUFFER),pos(NBUFFER,3),p(NBUFFER,3),stepl
real(8),allocatable,intent(in) :: atype(:),pos(:,:)
real(8),intent(in) :: p(mcx%NATOMS,3)
real(8),intent(in) :: stepl
real(8) :: potentialEnergy

!real(8) :: atypeTmp(NBUFFER),posTmp(NBUFFER,3),fTmp(NBUFFER,3),qTmp(NBUFFER)
! TODO: v is dummy but COPYATOM(MOVE) moves it. Thus needs to allocate 3xNBUFFER array here.
!real(8) :: vdummy(NBUFFER,3) 
real(8) :: GPE(0:13)
integer :: NATOMSTmp

avsTmp%pos(1:mcx%NATOMS,1:3)=pos(1:mcx%NATOMS,1:3)+stepl*p(1:mcx%NATOMS,1:3)
avsTmp%atype(1:mcx%NATOMS)=atype(1:mcx%NATOMS)

! TODO: local number of atoms will change after COPYATOM(MOVE). 
! Need to retrieve the origianl value consistent with the coordinates before the move.
NATOMSTmp = mcx%NATOMS

!call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0], atypeTmp, posTmp, vdummy, fTmp, qTmp)
!call QEq(atypeTmp, posTmp, qTmp)
!call FORCE(atypeTmp, posTmp, fTmp, qTmp)
call COPYATOMS(mcx, avsTmp, qvt_cg, mpt, MODE_MOVE, [0.d0, 0.d0, 0.d0])
call QEq(ffp, avsTmp, qvt_cg, mpt, rxp, mcx)
call FORCE(ffp, mpt, bos_cg, avsTmp, qvt_cg, mcx)

mcx%NATOMS = NATOMSTmp 

mcx%PE(0)=sum(mcx%PE(1:13))
call MPI_ALLREDUCE(mcx%PE, GPE, size(mcx%PE), MPI_DOUBLE_PRECISION, MPI_SUM, mpt%mycomm, mpt%ierr)
potentialEnergy=GPE(0)

end function EvaluateEnergyWithStep

!---------------------------------------------------------------------------------


end module  
