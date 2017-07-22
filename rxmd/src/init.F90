module init_funcs

contains

!------------------------------------------------------------------------------------------
SUBROUTINE INITSYSTEM(mcx, ffp, avs, qvt, bos, rxp, cla, mpt)
use md_context; use atom_vars; use rxmd_params; use cmdline_args; use mpi_vars
use ff_params; use bo; use energy_terms; use fileio_funcs; use MemoryAllocator
use support_funcs; use qeq_vars; use qeq_funcs
! This subroutine takes care of setting up initial system configuration.
! Unit conversion of parameters (energy, length & mass) are also done here.
!------------------------------------------------------------------------------------------
implicit none

type(md_context_type),intent(inout) :: mcx
type(mpi_var_type),intent(in) :: mpt
type(cmdline_arg_type),intent(in) :: cla
type(forcefield_params),intent(in) :: ffp 
type(atom_var_type),intent(inout) :: avs
type(qeq_var_type),intent(inout) :: qvt
type(bo_var_type),intent(inout) :: bos
type(rxmd_param_type),intent(inout) :: rxp

integer :: i,j,k, ity, l(3), ist=0, ierr
real(8) :: mm, gmm, dns, mat(3,3)
integer(8) :: i8
real(8) :: rcsize(3), maxrcell

call GetRxmdParams(rxp, cla%ParmPath)

!--- an error trap
if(rxp%vprocs(1)*rxp%vprocs(2)*rxp%vprocs(3) /= mpt%nprocs ) then
  if(mpt%myid==0) write(6,'(a60,3i3,i5)')  &
  "ERROR: requested/allocated # of procs are not consistent: ", rxp%vprocs(1:3), mpt%nprocs
  call MPI_FINALIZE(mpt%ierr)
  stop
endif

!--- initialize charge with QEq
if(rxp%mdmode==0) then
  if(mpt%myid==0) then
    print'(a,f12.3,a,i6,a)', &
         'INFO: mdmode==0, setting isQEQ is 1. Atomic velocities are scaled to ', &
          rxp%treq, ' [K] every ', rxp%sstep, ' steps.'
  endif
  rxp%isQEq=1
endif

!--- square the spring const in the extended Lagrangian method 
qvt%Lex_w2=2.d0*rxp%Lex_k/rxp%dt/rxp%dt

!--- setup the vector ID and parity for processes, in x, y and z order.
mcx%vID(1)=mod(mpt%myid,rxp%vprocs(1))
mcx%vID(2)=mod(mpt%myid/rxp%vprocs(1),rxp%vprocs(2))
mcx%vID(3)=mpt%myid/(rxp%vprocs(1)*rxp%vprocs(2))
mcx%myparity(1)=mod(mcx%vID(1),2)
mcx%myparity(2)=mod(mcx%vID(2),2)
mcx%myparity(3)=mod(mcx%vID(3),2)

k=0
do i=1,3
!--- Add (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1) to vector ID applying PBC. 
   do j=1,-1,-2

!--- save vector ID.
      l(1:3) = mcx%vID(1:3)

!--- apply PBC.
      l(i)=mod( (mcx%vID(i) + j + rxp%vprocs(i)), rxp%vprocs(i) )
     
!--- get direction index.
      k=k+1 ![123456] 

!--- convert vector ID to sequential ID.
      mcx%target_node(k) = l(1)+l(2)*rxp%vprocs(1)+l(3)*rxp%vprocs(1)*rxp%vprocs(2) 
   enddo         
                   
enddo    

!--- dt/2*mass, mass/2
call allocatord1d(mcx%dthm, 1, ffp%nso)
call allocatord1d(mcx%hmas, 1, ffp%nso)
do ity=1,ffp%nso
   mcx%dthm(ity) = rxp%dt*0.5d0/ffp%mass(ity)
   mcx%hmas(ity) = 0.5d0*ffp%mass(ity)
enddo

call initialize_atom_vars(avs, mcx%NBUFFER)
call ReadBIN(mcx, avs, qvt, rxp, mpt, trim(cla%dataDir)//"/rxff.bin")

call initialize_bo_vars(bos, mcx%NBUFFER)
call initialize_energy_terms(mcx%NBUFFER)
call initialize_qeq_vars(qvt, mcx%NBUFFER, MAXNEIGHBS10)

!--- get total number of atoms per type. This will be used to determine
!--- subroutine cutofflength() 
allocate(mcx%natoms_per_type(ffp%nso),mcx%ibuf8(ffp%nso)) ! NOTE 8byte int is not supported in MemoryAllocator
mcx%natoms_per_type(:)=0
do i=1, mcx%NATOMS
   ity=nint(avs%atype(i))
   mcx%natoms_per_type(ity)=mcx%natoms_per_type(ity)+1
enddo

call MPI_ALLREDUCE(mcx%natoms_per_type, mcx%ibuf8, ffp%nso, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)
mcx%natoms_per_type(1:ffp%nso)=mcx%ibuf8(1:ffp%nso)
deallocate(mcx%ibuf8)

!--- get global number of atoms
i8=mcx%NATOMS ! Convert 4 byte to 8 byte
call MPI_ALLREDUCE(i8, mcx%GNATOMS, 1, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- determine cutoff distances only for exsiting atom pairs
call CUTOFFLENGTH(mcx, ffp)

!--- update box-related variables
call UpdateBoxParams(mcx, rxp)

#ifdef STRESS
!--- stress variables
call allocatord2d(astr(1,6,1,mcx%NBUFFER)
astr(:,:)=0.d0; 
#endif

!--- Linked List & Near Neighb Parameters
call allocatori2d(mcx%nbrlist,1,mcx%NBUFFER,0,MAXNEIGHBS)
call allocatori2d(mcx%nbrindx,1,mcx%NBUFFER,1,MAXNEIGHBS)
call allocatori2d(mcx%nbplist,1,mcx%NBUFFER,0,MAXNEIGHBS10)
call allocatori1d(mcx%llist,1,mcx%NBUFFER)
call allocatori3d(mcx%header,-MAXLAYERS,mcx%cc(1)-1+MAXLAYERS,-MAXLAYERS,mcx%cc(2)-1+MAXLAYERS,-MAXLAYERS,mcx%cc(3)-1+MAXLAYERS)
call allocatori3d(mcx%nacell,-MAXLAYERS,mcx%cc(1)-1+MAXLAYERS,-MAXLAYERS,mcx%cc(2)-1+MAXLAYERS,-MAXLAYERS,mcx%cc(3)-1+MAXLAYERS)

!--- returning force index array 
call allocatori1d(mcx%frcindx,1,mcx%NBUFFER)

!--- setup potential table
call POTENTIALTABLE(mcx, ffp)

!--- get real size of linked list cell
rcsize(1) = mcx%lata/rxp%vprocs(1)/mcx%cc(1)
rcsize(2) = mcx%latb/rxp%vprocs(2)/mcx%cc(2)
rcsize(3) = mcx%latc/rxp%vprocs(3)/mcx%cc(3)
maxrcell = maxval(rcsize(1:3))

!--- setup 10[A] radius mesh to avoid visiting unecessary cells 
call GetNonbondingMesh(mcx, ffp, rxp)

!--- get density 
mm = 0.d0
do i=1, mcx%NATOMS
   ity = nint(avs%atype(i))
   mm = mm + ffp%mass(ity)
enddo
call MPI_ALLREDUCE (mm, gmm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
dns=gmm/mcx%MDBOX*UDENS

!--- allocate & initialize Array size Stat variables
i=rxp%ntime_step/rxp%pstep+1
allocate(mcx%maxas(i,nmaxas))
mcx%maxas(:,:)=0

!--- print out parameters and open data file
if(mpt%myid==0) then
   write(6,'(a)') "----------------------------------------------------------------"
   write(6,'(a30,i9,a3,i9)') "req/alloc # of procs:", rxp%vprocs(1)*rxp%vprocs(2)*rxp%vprocs(3), "  /",mpt%nprocs
   write(6,'(a30,3i9)')      "req proc arrengement:", rxp%vprocs(1),rxp%vprocs(2),rxp%vprocs(3)
   write(6,'(a30,a70)')      "parameter set:", mcx%FFDescript
   write(6,'(a30,es12.2)')   "time step[fs]:",rxp%dt*UTIME
   write(6,'(a30,i3, i10, i10)') "MDMODE CURRENTSTEP NTIMESTPE:", &
                                  rxp%mdmode, mcx%current_step, rxp%ntime_step
   write(6,'(a30,i6,es10.1,i6,i6)') "isQEq,QEq_tol,NMAXQEq,qstep:", &
                                     rxp%isQEq,rxp%QEq_tol,rxp%NMAXQEq,rxp%qstep
   write(6,'(a30,f8.3,f8.3)') 'Lex_fqs,Lex_k:',rxp%Lex_fqs,rxp%Lex_k
   write(6,'(a30,f12.3,f8.3,i9)') 'treq,vsfact,sstep:',rxp%treq*UTEMP0, rxp%vsfact, rxp%sstep
   write(6,'(a30,2i6)') 'fstep,pstep:', rxp%fstep,rxp%pstep
   write(6,'(a30,i24,i24)') "NATOMS GNATOMS:", mcx%NATOMS, mcx%GNATOMS
   write(6,'(a30,3f12.3)') "LBOX:",mcx%LBOX(1:3)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",mcx%HH(1:3,1,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",mcx%HH(1:3,2,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",mcx%HH(1:3,3,0)
   print'(a30,3f12.3)', 'lata,latb,latc:', mcx%lata,mcx%latb,mcx%latc
   print'(a30,3f12.3)', 'lalpha,lbeta,lgamma:', mcx%lalpha,mcx%lbeta,mcx%lgamma
   write(6,'(a30,3f10.4)') "density [g/cc]:",dns
   write(6,'(a30,3i6)')  '# of linkedlist cell:', mcx%cc(1:3)
   write(6,'(a30,f10.3,2x,3f10.2)') "maxrc, lcsize [A]:", mcx%maxrc, &
         mcx%lata/mcx%cc(1)/rxp%vprocs(1), &
         mcx%latb/mcx%cc(2)/rxp%vprocs(2), &
         mcx%latc/mcx%cc(3)/rxp%vprocs(3)
   write(6,'(a30,3i6)')  '# of linkedlist cell (NB):', mcx%nbcc(1:3)
   write(6,'(a30,3f10.2)') "lcsize [A] (NB):", &
         mcx%lata/mcx%nbcc(1)/rxp%vprocs(1), &
         mcx%latb/mcx%nbcc(2)/rxp%vprocs(2), &
         mcx%latc/mcx%nbcc(3)/rxp%vprocs(3)
   write(6,'(a30,2i6)') "MAXNEIGHBS, MAXNEIGHBS10:", MAXNEIGHBS,MAXNEIGHBS10
   write(6,'(a30,i6,i9)') "NMINCELL, NBUFFER:", NMINCELL, mcx%NBUFFER
   write(6,'(a30,3(a12,1x))') "FFPath, DataDir, ParmPath:", &
                          trim(cla%FFPath), trim(cla%dataDir), trim(cla%ParmPath)

   print'(a30 $)','# of atoms per type:'
   do ity=1,ffp%nso
      if(mcx%natoms_per_type(ity)>0) print'(i12,a2,i2 $)',mcx%natoms_per_type(ity),' -',ity
   enddo
   print*

   print'(a)', "----------------------------------------------------------------"
   write(6,'(a)')  &
   "nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)"

endif

END SUBROUTINE

!------------------------------------------------------------------------------------------
SUBROUTINE INITVELOCITY(mcx, ffp, atype, v, treq)
use ff_params; use md_context; use mpi_vars
! Generate gaussian distributed velocity as an initial value  using Box-Muller algorithm
!------------------------------------------------------------------------------------------
implicit none

real(8),intent(in) :: treq
type(md_context_type),intent(inout) :: mcx
type(forcefield_params),intent(in) :: ffp 
real(8) :: atype(mcx%NBUFFER)
real(8) :: v(3,mcx%NBUFFER)

integer :: i, k, ity, ierr
real(8) :: vv(2), vsqr, vsl, rndm(2)
real(8) :: vCM(3), GvCM(3), mm, Gmm
real(8) :: vfactor

!--- assign velocity to two atoms together with BM algoritm. 
!--- If <NATOMS> is odd, the <NATOMS> + 1 element will be the ignored in later calculations.

do i=1, mcx%NATOMS, 2

  do k=1,3 ! three directions
     !--- generate gaussian distributed velocity
     vsqr=0.d0
     do while ( (vsqr >= 1.d0) .or. (vsqr==0.d0) ) 
        call random_number(rndm)
        vv(1) = 2.d0 * rndm(1) - 1.d0
        vv(2) = 2.d0 * rndm(2) - 1.d0
        vsqr = vv(1)**2 + vv(2)**2
     enddo

     vsl = sqrt(-2.d0 * log(vsqr)/vsqr)
     v(i,k)   = vv(1)*vsl
     v(i+1,k) = vv(2)*vsl
  enddo
  
enddo

!--- get the local momentum and mass.
vCM(:)=0.d0;  mm = 0.d0
do i=1, mcx%NATOMS
   ity = nint(atype(i))
   vCM(1:3)=vCM(1:3) + ffp%mass(ity)*v(i,1:3)
   mm = mm + ffp%mass(ity)
enddo
 
call MPI_ALLREDUCE (vCM, GvCM, size(vCM), MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE (mm, Gmm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- get the global momentum
GvCM(:)=GvCM(:)/Gmm

!--- set the total momentum to be zero and get the current kinetic energy. 
mcx%KE = 0.d0
do i=1, mcx%NATOMS
   v(i,1:3) = v(i,1:3) - GvCM(1:3)

   ity = nint(atype(i))
   mcx%KE = mcx%KE + mcx%hmas(ity)*sum( v(i,1:3)*v(i,1:3) )
enddo

call MPI_ALLREDUCE (mcx%KE, mcx%GKE, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
mcx%GKE = mcx%GKE/mcx%GNATOMS

!--- scale the obtained velocity to get the initial temperature.
vfactor = sqrt(1.5d0*treq/mcx%GKE)
v(:,:) = vfactor * v(:,:)

end subroutine

!------------------------------------------------------------------------------------------
subroutine CUTOFFLENGTH(mcx, ffp)
use md_context; use ff_params 
!------------------------------------------------------------------------------------------
implicit none

type(md_context_type),intent(inout) :: mcx
type(forcefield_params),intent(in) :: ffp 

integer :: ity,jty,inxn
real(8) :: dr,BOsig

!--- get the cutoff length based on sigma bonding interaction.

! --- Remark --- 
! sigma bond before correction is the longest, namely longer than than pi and double-pi bonds
! thus check only sigma bond convergence.
allocate(mcx%rc(ffp%nboty), mcx%rc2(ffp%nboty))

mcx%rc(:)=0.d0; mcx%rc2(:)=0.d0
do ity=1,ffp%nso
do jty=ity,ffp%nso

   inxn=ffp%inxn2(ity,jty)
   if(inxn==0) cycle

   dr = 1.0d0
   BOsig=1.d0
   do while (BOsig > MINBOSIG) 
      dr = dr + 0.01d0
      BOsig = exp( ffp%pbo1(inxn)*(dr/ffp%r0s(ity,jty))**ffp%pbo2(inxn) ) !<- sigma bond prime
   enddo

   mcx%rc(inxn)  = dr
   mcx%rc2(inxn) = dr*dr
enddo
enddo

!----------------------------------------------------------------------------
! In some cases, an atom that do not exist in simulation gives
! the longest bond-order cutoff length. the check below is to ignore
! such atoms to keep the linkedlist cell dimensions as small as possible.
!----------------------------------------------------------------------------
do ity=1,ffp%nso
   if(mcx%natoms_per_type(ity)==0) then
      do jty=1,ffp%nso
         inxn=ffp%inxn2(ity,jty)
         if(inxn/=0) mcx%rc(inxn)=0.d0
         inxn=ffp%inxn2(jty,ity)
         if(inxn/=0) mcx%rc(inxn)=0.d0
      enddo
   endif
enddo

!--- get the max cutoff length 
mcx%maxrc=maxval(mcx%rc(:))

end subroutine

!------------------------------------------------------------------------------------------
subroutine POTENTIALTABLE(mcx, ffp)
use md_context; use ff_params; use MemoryAllocator
!------------------------------------------------------------------------------------------
implicit none

type(md_context_type),intent(inout) :: mcx
type(forcefield_params),intent(in) :: ffp 

integer :: i, ity,jty,inxn
real(8) :: dr1, dr2, dr3, dr4, dr5, dr6, dr7

real(8) :: exp1, exp2
real(8) :: gamwinvp, gamWij, alphaij, Dij0, rvdW0
real(8) :: Tap, dTap, fn13, dfn13, dr3gamij, rij_vd1

!--- first element in table 0: potential
!---                        1: derivative of potential
call allocatord3d(mcx%TBL_EClmb,0,1,1,NTABLE,1,ffp%nboty)
call allocatord3d(mcx%TBL_Evdw,0,1,1,NTABLE,1,ffp%nboty)
call allocatord2d(mcx%TBL_EClmb_QEq,1,NTABLE,1,ffp%nboty)

!--- unit distance in r^2 scale
mcx%UDR = ffp%rctap2/NTABLE
mcx%UDRi = 1.d0/mcx%UDR

do ity=1, ffp%nso
do jty=ity, ffp%nso

   inxn = ffp%inxn2(ity,jty)
   if(inxn/=0) then
      do i=1, NTABLE

         dr2 = mcx%UDR*i
         dr1 = sqrt(dr2)

!--- Interaction Parameters:
         gamWij = ffp%gamW(ity,jty)
         alphaij = ffp%alpij(ity,jty)
         Dij0 = ffp%Dij(ity,jty)
         rvdW0 = ffp%rvdW(ity,jty) 
         gamwinvp = (1.d0/gamWij)**ffp%pvdW1

         dr3 = dr1*dr2
         dr4 = dr2*dr2
         dr5 = dr1*dr2*dr2
         dr6 = dr2*dr2*dr2
         dr7 = dr1*dr2*dr2*dr2 

         rij_vd1 = dr2**ffp%pvdW1h
         Tap = ffp%CTap(7)*dr7 + ffp%CTap(6)*dr6 + &
               ffp%CTap(5)*dr5 + ffp%CTap(4)*dr4 + ffp%CTap(0)
         fn13 = (rij_vd1 + gamwinvp)**ffp%pvdW1inv
         exp1 = exp( alphaij*(1.d0 - fn13 / rvdW0) )
         exp2 = sqrt(exp1)

         dr3gamij = ( dr3 + ffp%gamij(ity,jty) )**( -1.d0/3.d0 )

!!--- Energy Calculation:
!      PEvd = Tap*Dij0*(exp1 - 2d0*exp2)      
!      PEclmb = Tap*Cclmb*q(i)*q(j)*dr3gamij
!if(myid==0) print*,i, Tap*Dij0*(exp1 - 2d0*exp2), Tap*Cclmb*dr3gamij

          mcx%TBL_Evdw(0,i,inxn) = Tap*Dij0*(exp1 - 2d0*exp2)      
          mcx%TBL_Eclmb(0,i,inxn) = Tap*Cclmb*dr3gamij
          mcx%TBL_Eclmb_QEq(i,inxn) = Tap*Cclmb0_qeq*dr3gamij

!if(inxn==1.and.myid==0) print*,i,TBL_Evdw(i,inxn,0), TBL_Eclmb(i,inxn,0)

!--- Force Calculation:
         dTap = 7d0*ffp%CTap(7)*dr5 + 6d0*ffp%CTap(6)*dr4 + &
                5d0*ffp%CTap(5)*dr3 + 4d0*ffp%CTap(4)*dr2

         dfn13 = ((rij_vd1 + gamwinvp)**(ffp%pvdW1inv-1.d0)) * (dr2**(ffp%pvdW1h-1.d0)) 

!      CEvdw = Dij0*( dTap*(exp1 - 2.d0*exp2)  &
!           - Tap*(alphaij/rvdW0)*(exp1 - exp2)*dfn13 )
!      CEclmb = Cclmb*q(i)*q(j)*dr3gamij*( dTap - (dr3gamij**3)*Tap*dr(0) ) 

         mcx%TBL_Evdw(1,i,inxn) = Dij0*( dTap*(exp1 - 2.d0*exp2)  &
                            - Tap*(alphaij/rvdW0)*(exp1 - exp2)*dfn13 )
         mcx%TBL_Eclmb(1,i,inxn) = Cclmb*dr3gamij*( dTap - (dr3gamij**3)*Tap*dr1 ) 

      enddo
   endif

enddo
enddo

end subroutine

!----------------------------------------------------------------
subroutine GetNonbondingMesh(mcx, ffp, rxp)
use md_context; use rxmd_params; use ff_params; use MemoryAllocator
! setup 10[A] radius mesh to avoid visiting unecessary cells 
!----------------------------------------------------------------
implicit none

type(md_context_type),intent(inout) :: mcx
type(forcefield_params),intent(in) :: ffp 
type(rxmd_param_type),intent(in) :: rxp

integer :: i,j,k

real(8) :: latticePerNode(3), rr(3), dr2
real(8) :: maxrcell
integer :: imesh(3), maximesh, ii(3), i1

!--- initial estimate of LL cell dims
mcx%nblcsize(1:3)=3d0

!--- get mesh resolution which is close to the initial value of rlc.
latticePerNode(1)=mcx%lata/rxp%vprocs(1)
latticePerNode(2)=mcx%latb/rxp%vprocs(2)
latticePerNode(3)=mcx%latc/rxp%vprocs(3)
mcx%nbcc(1:3)=int(latticePerNode(1:3)/mcx%nblcsize(1:3))
mcx%nblcsize(1:3)=latticePerNode(1:3)/mcx%nbcc(1:3)
maxrcell = maxval(mcx%nblcsize(1:3))

!--- get # of linked list cell to cover up the non-bonding (10[A]) cutoff length
imesh(1:3)  = int(ffp%rctap/mcx%nblcsize(1:3)) + 1
maximesh = maxval(imesh(1:3))

!--- List up only cell indices within the cutoff range.
!--- pre-compute nmesh to get exact array size.
mcx%nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   ii(1:3) = [i,j,k]
   do i1 = 1, 3
      if(ii(i1)>0) then
         ii(i1)=ii(i1)-1
      else if(ii(i1)<0) then
         ii(i1)=ii(i1)+1
      endif
   enddo
   rr(1:3) = ii(1:3)*mcx%nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= ffp%rctap**2) mcx%nbnmesh = mcx%nbnmesh + 1
enddo; enddo; enddo

call allocatori2d(mcx%nbmesh,1,3,1,mcx%nbnmesh)

mcx%nbmesh(:,:)=0
mcx%nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   ii(1:3) = [i,j,k]
   do i1 = 1, 3
      if(ii(i1)>0) then
         ii(i1)=ii(i1)-1
      else if(ii(i1)<0) then
         ii(i1)=ii(i1)+1
      endif
   enddo
   rr(1:3) = ii(1:3)*mcx%nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= ffp%rctap**2) then
      mcx%nbnmesh = mcx%nbnmesh + 1
      mcx%nbmesh(1:3,mcx%nbnmesh) = (/i, j, k/)
   endif
enddo; enddo; enddo

call allocatori1d(mcx%nbllist,1,mcx%NBUFFER)
call allocatori3d(mcx%nbheader, &
                -MAXLAYERS_NB,mcx%nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,mcx%nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,mcx%nbcc(3)-1+MAXLAYERS_NB)
call allocatori3d(mcx%nbnacell, &
                -MAXLAYERS_NB,mcx%nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,mcx%nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,mcx%nbcc(3)-1+MAXLAYERS_NB)

!--- normalize nblcsize, like lcsize.
mcx%nblcsize(1:3)=mcx%nblcsize(1:3)/(/mcx%lata,mcx%latb,mcx%latc/)

end subroutine

end module
