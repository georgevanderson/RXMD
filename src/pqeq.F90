!------------------------------------------------------------------------------
subroutine PQEq(atype, pos, q)
!use atoms
use pqeq_vars
use parameters
! Two vector electronegativity equilization routine
!
! The linkedlist cell size is determined by the cutoff length of bonding 
! interaction <rc> = 3A. Since the non-bonding interaction cutoff <Rcut> = 10A,
! need to take enough layers to calculate non-bonding interactoins.
!
!<Gnew>, <Gold> :: NEW and OLD squre norm of Gradient vector.
!<Est> :: ElectroSTatic energy
!-------------------------------------------------------------------------------
implicit none

real(8),intent(in) :: atype(NBUFFER), pos(NBUFFER,3)
real(8),intent(out) :: q(NBUFFER)

!real(8) :: fpqeq(NBUFFER)
real(8) :: vdummy(1,1), fdummy(1,1)

integer :: i,j,l2g
integer :: i1,j1,k1, nmax
real(8) :: Gnew(2), Gold(2) 
real(8) :: Est, GEst1, GEst2, g_h(2), h_hsh(2)
real(4) :: lmin(2)
real(8) :: buf(4), Gbuf(4)
real(8) :: ssum, tsum, mu
real(8) :: qsum, gqsum
real(8) :: QCopyDr(3)

call system_clock(i1,k1)

QCopyDr(1:3)=rctap/(/lata,latb,latc/)

!--- Initialize <s> vector with current charge and <t> vector with zero.
!--- isQEq==1 Normal QEq, isQEq==2 Extended Lagrangian method, DEFAULT skip QEq 
select case(isQEq)

!=== original QEq ===!
  case (1) 
!--- In the original QEq, fictitious charges are initialized with real charges
!--- and set zero.
    qsfp(1:NATOMS)=q(1:NATOMS)
    qsfv(1:NATOMS)=0.d0
!--- Initialization of the two vector QEq 
    qs(:)=0.d0
    qt(:)=0.d0
    qs(1:NATOMS)=q(1:NATOMS)
    nmax=NMAXQEq

!=== Extended Lagrangian method ===!
  case(2)
!--- charge mixing.
    qs(1:NATOMS)=Lex_fqs*qsfp(1:NATOMS)+(1.d0-Lex_fqs)*q(1:NATOMS)
!--- the same as the original QEq, set t vector zero
    qt(1:NATOMS)=0.d0
!--- just run one step
    nmax=1

!=== else, just return ===!
  case default
     return

end select

#ifdef QEQDUMP 
open(91,file="qeqdump"//trim(rankToString(myid))//".txt")
open(95,file="qdump"//trim(rankToString(myid))//".txt")
open(94,file="pqeqdump"//trim(rankToString(myid))//".txt")
#endif

#ifdef PQEQDUMP
open(92,file="pqeqdump_before"//trim(rankToString(myid))//".txt")
open(93,file="pqeqdump_after"//trim(rankToString(myid))//".txt")
#endif
!--- copy atomic coords and types from neighbors, used in qeq_initialize()
call COPYATOMS_SC(MODE_COPY_SC, QCopyDr, atype, pos, vdummy, fdummy, q)
call LINKEDLIST(atype, pos, nblcsize, nbheader, nbllist, nbnacell, nbcc, MAXLAYERS_NB)

call qeq_initialize()

#ifdef QEQDUMP 
do i=1, NATOMS
   write(95,'(i10,es25.15)') i,q(i)
   do j1=1,nbplist(i,0)
      j = nbplist(i,j1)
      !write(91,'(4i6,4es25.15)') -1, l2g(atype(i)),nint(atype(i)),l2g(atype(j)),hessian(j1,i)
      write(91,'(5i6,4es25.15)') nstep, i,j,l2g(atype(i)),l2g(atype(j)),hessian(j1,i)
   enddo
enddo
#endif

#ifdef PQEQDUMP 
do i=1, NATOMS + na/ne
   do j1=1,nbplist_sc(i,0)
      j = nbplist_sc(i,j1)
      !write(94,'(4i6,4es25.15)') -1, l2g(atype(i)),nint(atype(i)),l2g(atype(j)),hessian_sc(j1,i)
      !write(94,'(5i6,4es25.15)') nstep, i,j,l2g(atype(i)),l2g(atype(j)),hessian_sc(j1,i)
         if (l2g(atype(i)) < l2g(atype(j))) then
          write(94,'(5i6,4es25.15)') nstep, i,j,l2g(atype(i)),l2g(atype(j)),hessian_sc(j1,i)
       else
          write(94,'(5i6,4es25.15)') nstep, j,i,l2g(atype(j)),l2g(atype(i)),hessian_sc(j1,i)
       endif
   enddo
enddo
#endif

!--- after the initialization, only the normalized coords are necessary for COPYATOMS_SC()
!--- The atomic coords are converted back to real at the end of this function.
call COPYATOMS_SC(MODE_QCOPY1_SC,QCopyDr, atype, pos, vdummy, fdummy, q)
call get_gradient_sc(Gnew)

!--- Let the initial CG direction be the initial gradient direction
hs(1:NATOMS) = gs(1:NATOMS)
ht(1:NATOMS) = gt(1:NATOMS)

call COPYATOMS_SC(MODE_QCOPY2_SC,QCopyDr, atype, pos, vdummy, fdummy, q)

GEst2=1.d99
do nstep_qeq=0, nmax-1

#ifdef QEQDUMP 
  qsum = sum(q(1:NATOMS))
  call MPI_ALLREDUCE(qsum, gqsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
#endif

  !call get_hsh(Est)
  call get_hsh_sc(Est)

#ifdef PQEQDUMP 
do i=1, NATOMS
   !do j1=1,nbplist_sc(i,0)
   !   j = nbplist_sc(i,j1)
      write(92,'(i6,2es25.15)') i, hshs(i), hsht(i)
   !enddo
enddo

#endif

#ifdef DEBUG_CPBK
print '(a20,4es25.15)', "(hshs,hsht) before", sum(hshs(1:NATOMS)), sum(hsht(1:NATOMS)), &
      sum(hshs(NATOMS+1:NATOMS+na/ne)), sum(hsht(NATOMS+1:NATOMS+na/ne))
#endif

  call COPYATOMS_SC(MODE_CPHSH_SC,QCopyDr, atype, pos, vdummy, fdummy, q)

#ifdef DEBUG_CPBK
print '(a20,4es25.15)', "(hshs,hsht) after ", sum(hshs(1:NATOMS)), sum(hsht(1:NATOMS)), &
      sum(hshs(NATOMS+1:NATOMS+na/ne)), sum(hsht(NATOMS+1:NATOMS+na/ne))
#endif

  call MPI_ALLREDUCE(Est, GEst1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

#ifdef PQEQDUMP 
do i=1, NATOMS
   !do j1=1,nbplist_sc(i,0)
   !   j = nbplist_sc(i,j1)
      write(93,'(i6,2es25.15)') i, hshs(i), hsht(i)
   !enddo
enddo
#endif

#ifdef QEQDUMP
  if(myid==0) print'(i5,6es25.15)', nstep_qeq, 0.5d0*log(Gnew(1:2)/GNATOMS), GEst1, GEst2, GEst1-GEst2, gqsum
#endif
  if( ( 0.5d0*( abs(GEst2) + abs(GEst1) ) < QEq_tol) ) exit 
  if( abs(GEst2) > 0.d0 .and. (abs(GEst1/GEst2-1.d0) < QEq_tol) ) exit
  GEst2 = GEst1

!--- line minimization factor of <s> vector
  g_h(1) = dot_product(gs(1:NATOMS), hs(1:NATOMS))
  h_hsh(1) = dot_product(hs(1:NATOMS), hshs(1:NATOMS))

!--- line minimization factor of <t> vector
  g_h(2) = dot_product(gt(1:NATOMS), ht(1:NATOMS))
  h_hsh(2) = dot_product(ht(1:NATOMS), hsht(1:NATOMS))

  buf(1)=g_h(1);   buf(2)=g_h(2)
  buf(3)=h_hsh(1); buf(4)=h_hsh(2)
  Gbuf(:)=0.d0
  call MPI_ALLREDUCE(buf, Gbuf, 4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
  g_h(1) = Gbuf(1);   g_h(2) = Gbuf(2)
  h_hsh(1) = Gbuf(3); h_hsh(2) = Gbuf(4)

  lmin(1:2) = g_h(1:2)/h_hsh(1:2)

!--- line minimization for each vector
  qs(1:NATOMS) = qs(1:NATOMS) + lmin(1)*hs(1:NATOMS)
  qt(1:NATOMS) = qt(1:NATOMS) + lmin(2)*ht(1:NATOMS)

!--- get a current electronegativity <mu>
  ssum = sum(qs(1:NATOMS))
  tsum = sum(qt(1:NATOMS))
  buf(1) = ssum; buf(2) = tsum

  Gbuf(:)=0.d0
  call MPI_ALLREDUCE(buf, Gbuf, 2, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
  ssum=Gbuf(1); tsum=Gbuf(2)

  mu = ssum/tsum

!--- update atom charges
  q(1:NATOMS) = qs(1:NATOMS) - mu*qt(1:NATOMS)

!--- update new charges of buffered atoms.
  call COPYATOMS_SC(MODE_QCOPY1_SC,QCopyDr, atype, pos, vdummy, fdummy, q)

#ifdef DEBUG_CPBK
print '(a20,1es25.15)', "(q)",sum(q(1:NATOMS))
#endif

!--- save old residues.  
  Gold(:) = Gnew(:)
  call get_gradient_sc(Gnew)
  !call get_gradient(Gnew)

!--- get new conjugate direction
  hs(1:NATOMS) = gs(1:NATOMS) + (Gnew(1)/Gold(1))*hs(1:NATOMS)
  ht(1:NATOMS) = gt(1:NATOMS) + (Gnew(2)/Gold(2))*ht(1:NATOMS)

!--- update new conjugate direction for buffered atoms.
  call COPYATOMS_SC(MODE_QCOPY2_SC,QCopyDr, atype, pos, vdummy, fdummy, q)

enddo

!--- for PQEq
!call update_shell_positions()
call update_shell_positions_sc()

call system_clock(j1,k1)
it_timer(1)=it_timer(1)+(j1-i1)

! save # of QEq iteration 
it_timer(24)=it_timer(24)+nstep_qeq

#ifdef QEQDUMP 
close(91)
close(95)
#endif

#ifdef PQEQDUMP 
close(92)
close(93)
close(94)
#endif
return 

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------
subroutine update_shell_positions()
implicit none
!-----------------------------------------------------------------------------------------------------------------------
integer :: i,ity,j,jty,j1,inxn
real(8) :: shelli(3),shellj(3), qjc, clmb, dclmb, dr2
real(8) :: sforce(NATOMS,3), sf(3), Esc, Ess
real(8) :: ff(3)

real(8) :: dr(3)

sforce(1:NATOMS,1:3)=0.d0
do i=1, NATOMS

   ity = nint(atype(i))

   ! if i-atom is not polarizable, no force acting on i-shell. 
   if( .not. isPolarizable(ity) ) cycle 

   if(isEfield) sforce(i,eFieldDir) = sforce(i,eFieldDir) + Zpqeq(ity)*eFieldStrength*Eev_kcal

   sforce(i,1:3) = sforce(i,1:3) - Kspqeq(ity)*spos(i,1:3) ! Eq. (37)
   shelli(1:3) = pos(i,1:3) + spos(i,1:3)

   do j1 = 1, nbplist(i,0)

      j = nbplist(i,j1)
      jty = nint(atype(j))

      qjc = q(j) + Zpqeq(jty)
      shellj(1:3) = pos(j,1:3) + spos(j,1:3)

      ! j-atom can be either polarizable or non-polarizable. In either case,
      ! there will be force on i-shell from j-core.  qjc takes care of the difference.  Eq. (38)
      dr(1:3)=shelli(1:3)-pos(j,1:3)
      call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(ity,jty),Esc, inxnpqeq(ity, jty), TBL_Eclmb_psc,sf)

      ff(1:3)=-Cclmb0*sf(1:3)*qjc*Zpqeq(ity)
      sforce(i,1:3)=sforce(i,1:3)-ff(1:3)

      ! if j-atom is polarizable, there will be force on i-shell from j-shell. Eq. (38)
      if( isPolarizable(jty) ) then 
         dr(1:3)=shelli(1:3)-shellj(1:3)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphass(ity,jty),Ess, inxnpqeq(ity, jty), TBL_Eclmb_pss,sf)

         ff(1:3)=Cclmb0*sf(1:3)*Zpqeq(ity)*Zpqeq(jty)
         sforce(i,1:3)=sforce(i,1:3)-ff(1:3)

      endif

   enddo

enddo

!--- update shell positions after finishing the shell-force calculation.  Eq. (39)
do i=1, NATOMS
   ity = nint(atype(i))
   if( isPolarizable(ity) ) spos(i,1:3) = spos(i,1:3) + sforce(i,1:3)/Kspqeq(ity)
enddo


end subroutine

!-----------------------------------------------------------------------------------------------------------------------
subroutine update_shell_positions_sc()
implicit none
!-----------------------------------------------------------------------------------------------------------------------
integer :: i,ity,j,jty,j1,inxn
real(8) :: shelli(3),shellj(3), qic, qjc, clmb, dclmb, dr2
!real(8) :: sforce(NATOMS,3), sf(3), Esc, Ess
real(8) :: sf(3), Esci, Escj, Ess
real(8) :: ff(3)

real(8) :: dr(3)

sforce(1:NBUFFER,1:3)=0.d0

do i=1, NATOMS + na/ne

   ity = nint(atype(i))

   ! if i-atom is not polarizable, no force acting on i-shell. 
   !if( .not. isPolarizable(ity) ) cycle 


   !core-i shell-i interaction for resident atoms
   if( isPolarizable(ity) .and. i <= NATOMS ) then
      if(isEfield) sforce(i,eFieldDir) = sforce(i,eFieldDir) + Zpqeq(ity)*eFieldStrength*Eev_kcal
      sforce(i,1:3) = sforce(i,1:3) - Kspqeq(ity)*spos(i,1:3) ! Eq. (37)
   endif

   shelli(1:3) = pos(i,1:3) + spos(i,1:3)
   qic = q(i) + Zpqeq(ity)

   do j1 = 1, nbplist_sc(i,0)

      j = nbplist_sc(i,j1)
      jty = nint(atype(j))

      qjc = q(j) + Zpqeq(jty)
      shellj(1:3) = pos(j,1:3) + spos(j,1:3)

      ! j-atom can be either polarizable or non-polarizable. In either case,
      ! there will be force on i-shell from j-core.  qjc takes care of the difference.  Eq. (38)
      if( isPolarizable(ity) ) then
         dr(1:3)=shelli(1:3)-pos(j,1:3)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(ity,jty),Escj, inxnpqeq(ity, jty), TBL_Eclmb_psc,sf)

         ff(1:3)=-Cclmb0*sf(1:3)*qjc*Zpqeq(ity)
         sforce(i,1:3)=sforce(i,1:3)-ff(1:3)
      endif

      ! there will be force on j-shell from i-core.  qjc takes care of the difference.  Eq. (38)
      if( isPolarizable(jty) ) then
         dr(1:3)=shellj(1:3)-pos(i,1:3)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(jty,ity),Esci, inxnpqeq(jty, ity), TBL_Eclmb_psc,sf)

         ff(1:3)=-Cclmb0*sf(1:3)*qic*Zpqeq(jty)
         sforce(j,1:3)=sforce(j,1:3)-ff(1:3)
      endif

      ! if j-atom is polarizable, there will be force on i-shell from j-shell. Eq. (38)
      if(  isPolarizable(ity) .and. isPolarizable(jty) ) then 
         dr(1:3)=shelli(1:3)-shellj(1:3)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphass(ity,jty),Ess, inxnpqeq(ity, jty), TBL_Eclmb_pss,sf)

         ff(1:3)=Cclmb0*sf(1:3)*Zpqeq(ity)*Zpqeq(jty)
         sforce(i,1:3)=sforce(i,1:3)-ff(1:3)
         sforce(j,1:3)=sforce(j,1:3)+ff(1:3) !reaction force on shell-j
      endif

   enddo

enddo

!print *,"before shell update"
#ifdef DEBUG_CPBK
print '(a20,3es25.15)', "(sforce) before", sum(sforce(1:NATOMS,1)),sum(sforce(1:NATOMS,2)),sum(sforce(1:NATOMS,3))
#endif

call COPYATOMS_SC(MODE_CPBKSHELL_SC, QCopyDr, atype, pos, vdummy, fdummy, q)

#ifdef DEBUG_CPBK
print '(a20,3es25.15)', "(sforce) after", sum(sforce(1:NATOMS,1)),sum(sforce(1:NATOMS,2)),sum(sforce(1:NATOMS,3))
#endif
!print *,"after shell update"

!--- update shell positions after finishing the shell-force calculation.  Eq. (39)
do i=1, NATOMS
   ity = nint(atype(i))
   if( isPolarizable(ity) ) spos(i,1:3) = spos(i,1:3) + sforce(i,1:3)/Kspqeq(ity)
enddo


end subroutine


!-----------------------------------------------------------------------------------------------------------------------
subroutine qeq_initialize()
use atoms; use parameters; use MemoryAllocator
! This subroutine create a neighbor list with cutoff length = 10[A] and save the hessian into <hessian>.  
! <nbrlist> and <hessian> will be used for different purpose later.
!-----------------------------------------------------------------------------------------------------------------------
implicit none
integer :: i,j, ity, jty, n, m, mn, nn
integer :: c1,c2,c3, c4,c5,c6
integer :: ci1,ci2,ci3,ci4,ci5,ci6
real(4) :: dr2
real(8) :: dr(3), drtb
real(8) :: alphaij, pqeqc, pqeqs, ff(3)
integer :: itb, inxn

integer :: ti,tj,tk

logical :: sc_test,duplicate_test

call system_clock(ti,tk)

nbplist(:,0) = 0
nbplist_sc(:,0) = 0
fpqeq(:) = 0.d0

#if 0

!$omp parallel do schedule(runtime), default(shared), &
!$omp private(i,j,ity,jty,n,m,mn,nn,c1,c2,c3,c4,c5,c6,dr,dr2,drtb,itb,inxn,pqeqc,pqeqs,ff)
do c1=0, nbcc(1)-1
do c2=0, nbcc(2)-1
do c3=0, nbcc(3)-1

   i = nbheader(c1,c2,c3)
   do m = 1, nbnacell(c1,c2,c3)

   ity=nint(atype(i))

   !fpqeq(i)=0.d0

   do mn = 1, nbnmesh
      c4 = c1 + nbmesh(1,mn)
      c5 = c2 + nbmesh(2,mn)
      c6 = c3 + nbmesh(3,mn)


      j = nbheader(c4,c5,c6)
      do n=1, nbnacell(c4,c5,c6)

         if(i/=j) then
            dr(1:3) = pos(i,1:3) - pos(j,1:3)
            dr2 =  sum(dr(1:3)*dr(1:3))

            if(dr2 < rctap2) then

               jty = nint(atype(j))

!--- make neighbor-list upto the taper function cutoff
!$omp atomic
               nbplist(i,0) = nbplist(i,0) + 1
               nbplist(i,nbplist(i,0)) = j

!--- for SC, only one-way neighbor and neighbor atoms in + direction are needed
               !sc_test = (j > i .and. (pos(j,1) .ge. 0.d0) .and. (pos(j,2) .ge. 0.d0) .and. (pos(j,3) .ge. 0.d0)) 
               !if (sc_test .eqv. .TRUE.) then
               !   if(nbheader(c1,c2,c3) /= nbheader(c4,c5,c6) .or. (i<j)) then
!$omp atomic
               !   nbplist_sc(i,0) = nbplist_sc(i,0) + 1
               !   nbplist_sc(i,nbplist_sc(i,0)) = j
               !   endif
               !endif
!--- get table index and residual value
               itb = int(dr2*UDRi)
               drtb = dr2 - itb*UDR
               drtb = drtb*UDRi

!--- PEQq : 
               ! contribution from core(i)-core(j)
               call get_coulomb_and_dcoulomb_pqeq(dr,alphacc(ity,jty),pqeqc,inxnpqeq(ity, jty),TBL_Eclmb_pcc,ff)

               hessian(nbplist(i,0),i) = Cclmb0_qeq * pqeqc
               !if (sc_test .eqv. .TRUE.) then
               !   if(nbheader(c1,c2,c3) /= nbheader(c4,c5,c6) .or. (i<j)) then
               !      hessian_sc(nbplist_sc(i,0),i) = Cclmb0_qeq * pqeqc
               !   endif
               !endif

               fpqeq(i) = fpqeq(i) + Cclmb0_qeq * pqeqc * Zpqeq(jty) ! Eq. 30

               ! contribution from C(r_icjc) and C(r_icjs) if j-atom is polarizable
               if( isPolarizable(jty) ) then 
                  dr(1:3)=pos(i,1:3) - pos(j,1:3) - spos(j,1:3) ! pos(i,1:3)-(pos(j,1:3)+spos(j,1:3))  
                  call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(jty,ity),pqeqs,inxnpqeq(jty, ity),TBL_Eclmb_psc,ff)

                  fpqeq(i) = fpqeq(i) - Cclmb0_qeq * pqeqs * Zpqeq(jty) ! Eq. 30
               endif

            endif
         endif

         j=nbllist(j)
      enddo
   enddo !   do mn = 1, nbnmesh

   i=nbllist(i)
   enddo
enddo; enddo; enddo
!$omp end parallel do


#endif

!--- MATT TO DO-----------------------------
!--- Loop over cell near domain boundary to extract NT cell interactions (non-resident cell interaction)
!--- Determine NT interaction 
!print '(a,2i10)', "NATOMS,na/ne:",NATOMS, na/ne
!print '(a,3i10)', "nbcc(:):",nbcc(:)

!$omp parallel do schedule(runtime), default(shared), &
!$omp private(i,j,ity,jty,n,m,mn,nn,c1,c2,c3,ci1,ci2,ci3,ci4,ci5,ci6,dr,dr2,drtb,itb,inxn,pqeqc,pqeqs,ff)
do c1=0, nbcc(1)-1
do c2=0, nbcc(2)-1
do c3=0, nbcc(3)-1

do mn = 1, nbnmesh_sc

   ci1 = c1+nbmesh_sc(1,mn,1)
   ci2 = c2+nbmesh_sc(2,mn,1)
   ci3 = c3+nbmesh_sc(3,mn,1)
   ci4 = c1+nbmesh_sc(1,mn,2)
   ci5 = c2+nbmesh_sc(2,mn,2)
   ci6 = c3+nbmesh_sc(3,mn,2)

!print '(a,9i5)',"c1,c2,c3,ci1,ci2,ci3,ci4,ci5,ci6:",c1,c2,c3,ci1,ci2,ci3,ci4,ci5,ci6
!--- check if cell (ci1,ci2,ci3) or (ci4,ci5,ci6) are resident cells. 
!--- If at least one of them is resident cell, skipped since it is already computed in the previous loop
!   if ((ci1 < nbcc(1) .and. ci2 < nbcc(2) .and. ci3 < nbcc(3)) .or. &
!       (ci4 < nbcc(1) .and. ci5 < nbcc(2) .and. ci6 < nbcc(3))) then
!#ifdef MATT_DEBUG
!          print '(a,i4,i4,i4,i4,i4,i4)',"Cycle:",ci1,ci2,ci3,ci4,ci5,ci6
!#endif
!          cycle !--one of the interacting cells are resident cells. Skipped
!   endif

!#if MATT_DEBUG > 1
!          print '(a,i4,i4,i4,i4,i4,i4)',"Not Cycle:",ci1,ci2,ci3,ci4,ci5,ci6
!#endif
      
   i = nbheader(ci1,ci2,ci3)
   do m = 1, nbnacell(ci1,ci2,ci3)

   ity=nint(atype(i))

   !fpqeq(i)=0.d0

      j = nbheader(ci4,ci5,ci6)
      do n=1, nbnacell(ci4,ci5,ci6)

            !if ((j > NATOMS)) then 
            !   print '(a,2i6)',"(i,j)J>NATOMS:",i,j
            !if (i > NATOMS) then 
            !   print '(a,2i6)',"(i,j)I,J>NATOMS:",i,j
            !endif
            !endif

         !if(i/=j) then
         !compute only atoms in the same cell
         if((nbheader(ci1,ci2,ci3) /= nbheader(ci4,ci5,ci6)) .or. (i<j)) then
            dr(1:3) = pos(i,1:3) - pos(j,1:3)
            dr2 =  sum(dr(1:3)*dr(1:3))

            !print '(a,2i6,2es15.5)',"(i,j,dist,rctap2):",i,j,dr2,rctap2
            
            if(dr2 < rctap2) then

               jty = nint(atype(j))

!--- for SC, only one-way neighbor and neighbor atoms in + direction are needed
!$omp atomic
               nbplist_sc(i,0) = nbplist_sc(i,0) + 1
               nbplist_sc(i,nbplist_sc(i,0)) = j
!--- get table index and residual value
            itb = int(dr2*UDRi)
            drtb = dr2 - itb*UDR
            drtb = drtb*UDRi

!--- PEQq : 
            ! contribution from core(i)-core(j)
            call get_coulomb_and_dcoulomb_pqeq(dr,alphacc(ity,jty),pqeqc,inxnpqeq(ity, jty),TBL_Eclmb_pcc,ff)
            hessian_sc(nbplist_sc(i,0),i) = Cclmb0_qeq * pqeqc


            fpqeq(i) = fpqeq(i) + Cclmb0_qeq * pqeqc * Zpqeq(jty) ! Eq. 30
            fpqeq(j) = fpqeq(j) + Cclmb0_qeq * pqeqc * Zpqeq(ity) ! Eq. 30

            ! contribution from C(r_icjc) and C(r_icjs) if j-atom is polarizable
            if( isPolarizable(jty) ) then 
               dr(1:3)=pos(i,1:3) - pos(j,1:3) - spos(j,1:3) ! pos(i,1:3)-(pos(j,1:3)+spos(j,1:3))  
               call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(jty,ity),pqeqs,inxnpqeq(jty, ity),TBL_Eclmb_psc,ff)

               fpqeq(i) = fpqeq(i) - Cclmb0_qeq * pqeqs * Zpqeq(jty) ! Eq. 30
            endif

            if( isPolarizable(ity) ) then 
               dr(1:3)=pos(j,1:3) - pos(i,1:3) - spos(i,1:3) ! pos(i,1:3)-(pos(j,1:3)+spos(j,1:3))  
               call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(ity,jty),pqeqs,inxnpqeq(ity, jty),TBL_Eclmb_psc,ff)

               fpqeq(j) = fpqeq(j) - Cclmb0_qeq * pqeqs * Zpqeq(ity) ! Eq. 30
            endif

          endif !dr2 < rctap2
          
          endif !(nbheader(ci1,ci2,ci3) /= nbheader(ci4,ci5,ci6)) .or. (i<j)

         j=nbllist(j)
         enddo !loop atom in cell j
      i=nbllist(i)
      enddo !loop atom in cell i
   
   enddo !do mn = nbmesh_sc
enddo; enddo; enddo
!$omp end parallel do


#ifdef DEBUG_CPBK
print '(a20,2es25.15)', "(fpqeq) before", sum(fpqeq(1:NATOMS)),sum(fpqeq(NATOMS+1:NATOMS+na/ne))
#endif
call COPYATOMS_SC(MODE_CPFPQEQ_SC, QCopyDr, atype, pos, vdummy, fdummy, q)

#ifdef DEBUG_CPBK
print '(a20,2es25.15)', "(fpqeq) after ", sum(fpqeq(1:NATOMS)),sum(fpqeq(NATOMS+1:NATOMS+na/ne))
#endif

!-------END MATT----------------------------

!--- for array size stat
if(mod(nstep,pstep)==0) then
  nn=maxval(nbplist(1:NATOMS,0))
  i=nstep/pstep+1
  maxas(i,3)=nn
endif


#ifdef MATT_DEBUG
!--- Compare nbplist and nbplist_sc stat
!do i = 1,NBUFFER
!   do j = 1,nbplist(i,0)
!      print '(a,i10,i10)', "nbplist: ", i,nbplist(i,j)
!   enddo
!enddo

!do i = 1,NBUFFER
!   do j = 1,nbplist_sc(i,0)
!      print '(a,i10,i10)', "nbplist_sc: ", i,nbplist_sc(i,j)
!   enddo
!enddo

print *,"*******************************************"
print '(a,i10)', "Total nbplist",sum(nbplist(:,0))
print '(a,i10)', "Total nbplist_sc",sum(nbplist_sc(:,0))
print *,"*******************************************"
#endif


call system_clock(tj,tk)
it_timer(16)=it_timer(16)+(tj-ti)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_hsh_sc(Est)
use atoms; use parameters
! This subroutine updates hessian*cg array <hsh> and the electrostatic energy <Est>. The computing is based on SC algorithm  
!-----------------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(OUT) :: Est
integer :: i,j,j1, ity, jty, inxn
real(8) :: eta_ity, Est1, dr2, dr(3)

real(8) :: Ccicj,Csicj,Csjci,Csisj,shelli(3),shellj(3),qic,qjc,ff(3)
real(8) :: Eshell

integer :: ti,tj,tk
call system_clock(ti,tk)

hshs(:) = 0.d0
hsht(:) = 0.d0

Est = 0.d0
!$omp parallel do default(shared), reduction(+:Est) &
!$omp private(i,j,j1,ity,jty,eta_ity,Est1,Eshell,Ccicj,Csicj,Csisj,shelli,shellj,qic,qjc,ff,dr,dr2)
do i=1, NATOMS + na/ne
   ity = nint(atype(i))
   eta_ity = eta(ity)

   !avoid hshs and hsht double counting by compute only for resident atoms
   if (i <= NATOMS) then
      hshs(i) = hshs(i) + eta_ity*hs(i)
      hsht(i) = hsht(i) + eta_ity*ht(i)
   endif

!--- for PQEq
   qic = q(i) + Zpqeq(ity)
   shelli(1:3) = pos(i,1:3) + spos(i,1:3)

   !can be removed?
   !dr2 = sum(spos(i,1:3)*spos(i,1:3)) ! distance between core-and-shell for i-atom

   !avoid Est q(i)*q(i) double counting. Only do this for resident atoms
   if (i <= NATOMS) then

      !Eshell = 0.d0
      !if(isPolarizable(ity)) then
      !   dr2 = sum( spos(i,1:3)*spos(i,1:3) )
      !   Eshell = 0.5d0*Kspqeq(ity)*dr2
      !endif

      !Est = Est + CEchrge*(chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i)) + Eshell
      !Est = Est + CEchrge*(chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i))
      !Est = Est + CEchage*(chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i))
      Est = Est + chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i)
   endif

   do j1 = 1, nbplist_sc(i,0)
      j = nbplist_sc(i,j1)
      jty = nint(atype(j))

!--- for PQEq
      qjc = q(j) + Zpqeq(jty)
      shellj(1:3) = pos(j,1:3) + spos(j,1:3)

      Ccicj = 0.d0; Csicj=0.d0; Csisj=0.d0; Csjci=0.d0

!--- core-i/core-j
      !Ccicj = hessian_sc(j1,i)*qic*qjc
      !Ccicj = Ccicj*qic*qjc  !core i-j full
      Ccicj = hessian_sc(j1,i)*qic*qjc  !core i-j full
      !Ccicj = Ccicj*qic*qjc*0.5d0  !core i-j full
      !dr(1:3)=pos(i,1:3)-pos(j,1:3)
      !call get_coulomb_and_dcoulomb_pqeq(dr,alphacc(ity,jty),Ccicj,inxnpqeq(ity,jty),TBL_Eclmb_pcc,ff)
      !Ccicj = 0.5d0*Cclmb0_qeq*Ccicj*qic*qjc  !core i-j full
      !Ccicj = Cclmb0_qeq*Ccicj*qic*qjc  !core i-j full
      !Ccicj = Ccicj*Zpqeq(ity)*Zpqeq(jty)  !core i-j full
      !print '(a10,es25.15)',"Ccicj:",Ccicj

      !--- shell-i/core-j
      if(isPolarizable(ity)) then
         dr(1:3)=shelli(1:3)-pos(j,1:3)
         !print '(a10,6es25.15)',"dist:",shelli(:),pos(j,:)
         !print '(a10,3es25.15)',"dr:",dr(:)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(ity,jty),Csicj,inxnpqeq(ity,jty),TBL_Eclmb_psc,ff)
         !Csicj=-Cclmb0_qeq*Csicj*qic*Zpqeq(jty)
         Csicj=-Cclmb0_qeq*Csicj*qjc*Zpqeq(ity)
      endif
       
      !--- core-i/shell-j
      if(isPolarizable(jty)) then
         dr(1:3)=shellj(1:3)-pos(i,1:3)
         !print '(a10,6es25.15)',"dist:",shellj(:),pos(i,:)
         !print '(a10,3es25.15)',"dr:",dr(:)
         !call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(jty,ity),Csjci,inxnpqeq(jty,ity),TBL_Eclmb_psc,ff)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(jty,ity),Csjci,inxnpqeq(jty,ity),TBL_Eclmb_psc,ff)
         !Csjci=-Cclmb0_qeq*Csjci*qjc*Zpqeq(ity)
         Csjci=-Cclmb0_qeq*Csjci*qic*Zpqeq(jty)
      endif

         !--- shell-i/shell-j
      if(isPolarizable(ity) .and. isPolarizable(jty)) then
         dr(1:3)=shelli(1:3)-shellj(1:3)
         !print '(a10,6es25.15)',"dist:",shelli(:),shellj(:)
         !print '(a10,3es25.15)',"dr:",dr(:)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphass(ity,jty),Csisj,inxnpqeq(ity,jty),TBL_Eclmb_pss,ff)
         !Csisj=Cclmb0_qeq*Csisj*Zpqeq(ity)*Zpqeq(jty)
         Csisj=Cclmb0_qeq*Csisj*Zpqeq(ity)*Zpqeq(jty)
         !print '(a10,es25.15)',"Csisj:",Csisj
      endif

      hshs(i) = hshs(i) + hessian_sc(j1,i)*hs(j)
      hsht(i) = hsht(i) + hessian_sc(j1,i)*ht(j)

      hshs(j) = hshs(j) + hessian_sc(j1,i)*hs(i)
      hsht(j) = hsht(j) + hessian_sc(j1,i)*ht(i)


!--- In FS PQeQ, get half of potential energy, then sum it up if atoms are resident.
!--- But in SC PQeQ, get full potential energy (no duplicate atom pair) and sum up all atoms including cached atoms.
      !Est1 = 2.d0*Ccicj + Csicj + Csjci + 2.d0*Csisj
      Est1 = Ccicj + Csicj + Csjci + Csisj
      
      !if (i < j) then  
         !print '(a,2i5,es25.15)',"Energy (i,j):",i,j,Est1
         !print '(a,2i5,5es15.7)',"i,j,Ccicj,Csicj,Csjci,Csisj:",i,j,Est1,Ccicj,Csicj,Csjci,Csisj
      !else
         !print '(a,2i5,es25.15)',"Energy (i,j):",j,i,Est1
         !print '(a,2i5,5es15.7)',"i,j,Ccicj,Csicj,Csjci,Csisj:",j,i,Est1,Ccicj,Csicj,Csjci,Csisj
      !endif 

      !Est = Est + 0.5d0*Est1
      !if(j<=NATOMS) Est = Est + 0.5d0*Est1
      Est = Est + Est1
     
      !print '(a,i5,i5,r*4)',"get_hsh_sc(i,j),Ccicj,Csicj,Csjci,Csisj:",i,j
      !print *,"get_hsh_sc(i,j),Ccicj,Csicj,Csjci,Csisj:",i,j,Ccicj,Csicj,Csjci,Csisj
   enddo

enddo
!$omp end parallel do

call system_clock(tj,tk)
it_timer(18)=it_timer(18)+(tj-ti)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_hsh(Est)
use atoms; use parameters
! This subroutine updates hessian*cg array <hsh> and the electrostatic energy <Est>.  
!-----------------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(OUT) :: Est
integer :: i,j,j1, ity, jty, inxn
real(8) :: eta_ity, Est1, dr2, dr(3)

real(8) :: Ccicj,Csicj,Csisj,shelli(3),shellj(3),qic,qjc,ff(3)
real(8) :: Eshell

integer :: ti,tj,tk
call system_clock(ti,tk)

Est = 0.d0
!$omp parallel do default(shared), reduction(+:Est) &
!$omp private(i,j,j1,ity,jty,eta_ity,Est1,Eshell,Ccicj,Csicj,Csisj,shelli,shellj,qic,qjc,ff,dr,dr2)
do i=1, NATOMS
   ity = nint(atype(i))
   eta_ity = eta(ity)

   hshs(i) = eta_ity*hs(i)
   hsht(i) = eta_ity*ht(i)

!--- for PQEq
   qic = q(i) + Zpqeq(ity)
   shelli(1:3) = pos(i,1:3) + spos(i,1:3)

   dr2 = sum(spos(i,1:3)*spos(i,1:3)) ! distance between core-and-shell for i-atom

   Est = Est + chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i)

   do j1 = 1, nbplist(i,0)
      j = nbplist(i,j1)
      jty = nint(atype(j))

!--- for PQEq
      qjc = q(j) + Zpqeq(jty)
      shellj(1:3) = pos(j,1:3) + spos(j,1:3)

      Ccicj = 0.d0; Csicj=0.d0; Csisj=0.d0

      Ccicj = hessian(j1,i)
      Ccicj = Ccicj*qic*qjc*0.5d0

      if(isPolarizable(ity)) then
         dr(1:3)=shelli(1:3)-pos(j,1:3)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(ity,jty),Csicj,inxnpqeq(ity,jty),TBL_Eclmb_psc,ff)
         Csicj=-Cclmb0_qeq*Csicj*qic*Zpqeq(jty)

         if(isPolarizable(jty)) then
             dr(1:3)=shelli(1:3)-shellj(1:3)
             call get_coulomb_and_dcoulomb_pqeq(dr,alphass(ity,jty),Csisj,inxnpqeq(ity,jty),TBL_Eclmb_pss,ff)
             Csisj=Cclmb0_qeq*Csisj*Zpqeq(ity)*Zpqeq(jty)
         endif
      endif

      hshs(i) = hshs(i) + hessian(j1,i)*hs(j)
      hsht(i) = hsht(i) + hessian(j1,i)*ht(j)

!--- get half of potential energy, then sum it up if atoms are resident.
      Est1 = Ccicj + Csicj + Csisj

      Est = Est + 0.5d0*Est1
      if(j<=NATOMS) Est = Est + 0.5d0*Est1
   enddo

enddo
!$omp end parallel do

call system_clock(tj,tk)
it_timer(18)=it_timer(18)+(tj-ti)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_gradient_sc(Gnew)
use atoms; use parameters
! Update gradient vector <g> and new residue <Gnew>
!-----------------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(OUT) :: Gnew(2)
real(8) :: eta_ity, ggnew(2)
integer :: i,j,j1, ity

!real(8) :: gssum, gtsum
integer :: ti,tj,tk
call system_clock(ti,tk)

gssum(:) = 0.d0
gtsum(:) = 0.d0

!gs(:) = 0.d0
!gt(:) = 0.d0

!$omp parallel do default(shared), schedule(runtime), private(gssum, gtsum, eta_ity,i,j,j1,ity)
do i=1,NATOMS + na/ne

   do j1=1, nbplist_sc(i,0) 
      j = nbplist_sc(i,j1)
      gssum(i) = gssum(i) + hessian_sc(j1,i)*qs(j)
      gtsum(i) = gtsum(i) + hessian_sc(j1,i)*qt(j)

      gssum(j) = gssum(j) + hessian_sc(j1,i)*qs(i)
      gtsum(j) = gtsum(j) + hessian_sc(j1,i)*qt(i)
   enddo

enddo 
!$omp end parallel do

#ifdef DEBUG_CPBK
print '(a20,2es25.15)', "(gssum,gtsum) before", sum(gssum(1:NATOMS)),sum(gtsum(1:NATOMS))
#endif

call COPYATOMS_SC(MODE_CPGSGT_SC, QCopyDr, atype, pos, vdummy, fdummy, q)

#ifdef DEBUG_CPBK
print '(a20,2es25.15)', "(gssum,gtsum) after ", sum(gssum(1:NATOMS)),sum(gtsum(1:NATOMS))
#endif

do i=1,NATOMS

   ity = nint(atype(i))
   eta_ity = eta(ity)

   gs(i) = - chi(ity) - eta_ity*qs(i) - gssum(i) - fpqeq(i)
   gt(i) = - 1.d0     - eta_ity*qt(i) - gtsum(i)

enddo

ggnew(1) = dot_product(gs(1:NATOMS), gs(1:NATOMS))
ggnew(2) = dot_product(gt(1:NATOMS), gt(1:NATOMS))
call MPI_ALLREDUCE(ggnew, Gnew, size(ggnew), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)


call system_clock(tj,tk)
it_timer(19)=it_timer(19)+(tj-ti)


end subroutine


!-----------------------------------------------------------------------------------------------------------------------
subroutine get_gradient(Gnew)
use atoms; use parameters
! Update gradient vector <g> and new residue <Gnew>
!-----------------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(OUT) :: Gnew(2)
real(8) :: eta_ity, ggnew(2)
integer :: i,j,j1, ity

real(8) :: lgssum, lgtsum

integer :: ti,tj,tk
call system_clock(ti,tk)

!$omp parallel do default(shared), schedule(runtime), private(gssum, gtsum, eta_ity,i,j,j1,ity)
do i=1,NATOMS

   lgssum=0.d0
   lgtsum=0.d0
   do j1=1, nbplist(i,0) 
      j = nbplist(i,j1)
      lgssum = lgssum + hessian(j1,i)*qs(j)
      lgtsum = lgtsum + hessian(j1,i)*qt(j)
   enddo

   ity = nint(atype(i))
   eta_ity = eta(ity)

   gs(i) = - chi(ity) - eta_ity*qs(i) - lgssum - fpqeq(i)
   gt(i) = - 1.d0     - eta_ity*qt(i) - lgtsum

enddo 
!$omp end parallel do

ggnew(1) = dot_product(gs(1:NATOMS), gs(1:NATOMS))
ggnew(2) = dot_product(gt(1:NATOMS), gt(1:NATOMS))
call MPI_ALLREDUCE(ggnew, Gnew, size(ggnew), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

call system_clock(tj,tk)
it_timer(19)=it_timer(19)+(tj-ti)

end subroutine

end subroutine PQEq
