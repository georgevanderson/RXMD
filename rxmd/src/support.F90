module support_funcs

contains

!------------------------------------------------------------------------------
subroutine vkick(mcx, dtf, atype, v, f)
use md_context
!------------------------------------------------------------------------------
implicit none

type(md_context_type),intent(in) :: mcx
real(8),allocatable :: atype(:),v(:,:),f(:,:)

integer :: i, ity
real(8) :: dtf

do i=1,mcx%NATOMS
   ity = nint(atype(i))
   v(i,1:3) = v(i,1:3) + dtf*mcx%dthm(ity)*f(i,1:3)
enddo

end subroutine

!----------------------------------------------------------------------------------------
subroutine PRINTE(mcx, mpt, pstep, atype, v, q)
use md_context; use ff_params; use MemoryAllocator; use mpi_vars
! calculate the kinetic energy and sum up all of potential energies, then print them.
!----------------------------------------------------------------------------------------
implicit none

type(md_context_type) :: mcx
type(mpi_var_type) :: mpt
integer,intent(in) :: pstep
real(8),intent(in),allocatable :: atype(:), q(:), v(:,:)

integer :: i,ity,cstep, ierr
real(8),save :: wt0
real(8) :: qq=0.d0,tt=0.d0,ss=0.d0,buf(0:20),Gbuf(0:20)

i=mcx%nstep/pstep+1
mcx%maxas(i,1)=mcx%NATOMS

mcx%KE=0.d0
do i=1, mcx%NATOMS
   ity=nint(atype(i))
   mcx%KE = mcx%KE + mcx%hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo
qq=sum(q(1:mcx%NATOMS))

#ifdef STRESS
!--- pressure 
ss=sum(astr(1:3,1:mcx%NATOMS))
#endif

!--- potential energy 
mcx%PE(0)=sum(mcx%PE(1:13))

!--- copy data into buffer
buf(0:13) = mcx%PE(0:13)
buf(14) = mcx%KE; buf(15) = ss; buf(16) = qq
call MPI_ALLREDUCE (buf, Gbuf, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, mpt%mycomm, ierr)

!--- copy data from buffer
mcx%GPE(0:13) = Gbuf(0:13)
mcx%GKE = Gbuf(14); ss = Gbuf(15); qq = Gbuf(16)

!--- compute properties
mcx%GPE(:)=mcx%GPE(:)/mcx%GNATOMS
mcx%GKE=mcx%GKE/mcx%GNATOMS
tt=mcx%GKE*UTEMP
#ifdef STRESS
ss=ss/3.d0/MDBOX*USTRS
#endif 

!--- total energy
mcx%GTE = mcx%GKE + mcx%GPE(0)
if(mpt%myid==0) then
   
   cstep = mcx%nstep + mcx%current_step 

   write(6,'(i9,3es13.5,6es11.3,1x,3f8.2,i4,f8.2,f8.2)') cstep,mcx%GTE,mcx%GPE(0),mcx%GKE, &
   mcx%GPE(1),sum(mcx%GPE(2:4)),sum(mcx%GPE(5:7)),sum(mcx%GPE(8:9)),mcx%GPE(10),sum(mcx%GPE(11:13)), &
   tt, ss, qq, mcx%nstep_qeq, GetTotalMemory()*1e-9, MPI_WTIME()-wt0 

#ifdef STRESS
   write(6,'(6es13.5)') pint(1,1)*USTRS, pint(2,2)*USTRS, pint(3,3)*USTRS, &
                        pint(2,3)*USTRS, pint(3,1)*USTRS, pint(1,2)*USTRS
#endif

endif

!--- save current time
wt0 = MPI_WTIME()
end subroutine
!----------------------------------------------------------------------
subroutine angular_momentum(mcx, mass, atype, pos, v)
use md_context; use ff_params; use mpi_vars
!----------------------------------------------------------------------
implicit none

type(md_context_type),intent(in) :: mcx
real(8),allocatable,dimension(:) :: mass
real(8),allocatable :: atype(:), pos(:,:),v(:,:)

integer :: i,ity,ierr
real(8) :: com(3), Gcom(3), intsr(3,3), Gintsr(3,3), intsr_i(3,3), angm(3), Gangm(3), angv(3), mm, Gmm
real(8) :: dr(3), dv(3)

!--- get center of mass
com(:)=0.d0;     Gcom(:)=0.d0
mm=0.d0; Gmm=0.d0

do i=1, mcx%NATOMS
   ity = nint(atype(i))
   mm = mm + mass(ity)
   com(1:3) = mass(ity)*pos(i,1:3)
enddo

call MPI_ALLREDUCE(mm, Gmm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(com, Gcom, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gcom(1:3) = Gcom(1:3)/Gmm


!--- get the angular momentum and inertia tensor from the com
angm(:)=0.d0;    Gangm(:)=0.d0
intsr(:,:)=0.d0; Gintsr(:,:)=0.d0

do i=1, mcx%NATOMS
   dr(1:3) = pos(i,1:3) - Gcom(1:3)
   
   angm(1) = mass(ity)*( dr(2)*v(i,3)-dr(3)*v(i,2) )
   angm(2) = mass(ity)*( dr(3)*v(i,1)-dr(1)*v(i,3) )
   angm(3) = mass(ity)*( dr(1)*v(i,2)-dr(2)*v(i,1) )

   intsr(1,1) = mass(ity)*( dr(2)**2+dr(3)**2 )
   intsr(2,2) = mass(ity)*( dr(3)**2+dr(1)**2 )
   intsr(3,3) = mass(ity)*( dr(1)**2+dr(2)**2 )

   intsr(1,2) =-mass(ity)*( dr(1)*dr(2) )
   intsr(1,3) =-mass(ity)*( dr(1)*dr(3) )
   intsr(2,3) =-mass(ity)*( dr(2)*dr(3) )

   intsr(2,1) = intsr(1,2)
   intsr(3,1) = intsr(1,3)
   intsr(3,2) = intsr(2,3)
enddo

call MPI_ALLREDUCE(angm, Gangm, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(intsr, Gintsr, 9, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- get angular velocity
call matinv(Gintsr, intsr_i)

angv(1) = sum(intsr_i(1,1:3)*angm(1:3))
angv(2) = sum(intsr_i(2,1:3)*angm(1:3))
angv(3) = sum(intsr_i(3,1:3)*angm(1:3))


!--- correct rotational motion wrt CoM.
do i=1,mcx%NATOMS
   dr(1:3) = pos(i,1:3) - Gcom(1:3)
   dv(1) = angv(2)*dr(3) - angv(3)*dr(2)
   dv(2) = angv(3)*dr(1) - angv(1)*dr(3)
   dv(3) = angv(1)*dr(2) - angv(2)*dr(1)

   v(i,1:3) = v(i,1:3) - dv(1:3)
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine matinv(m1,m2)
! get inverse of m1 and save to m2
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8) :: m1(3,3), m2(3,3), detm

m2(1,1) = m1(2,2)*m1(3,3)-m1(2,3)*m1(3,2)
m2(1,2) = m1(1,3)*m1(3,2)-m1(1,2)*m1(3,3)
m2(1,3) = m1(1,2)*m1(2,3)-m1(1,3)*m1(2,2)
m2(2,1) = m1(2,3)*m1(3,1)-m1(2,1)*m1(3,3)
m2(2,2) = m1(1,1)*m1(3,3)-m1(1,3)*m1(3,1)
m2(2,3) = m1(1,3)*m1(2,1)-m1(1,1)*m1(2,3)
m2(3,1) = m1(2,1)*m1(3,2)-m1(2,2)*m1(3,1)
m2(3,2) = m1(1,2)*m1(3,1)-m1(1,1)*m1(3,2)
m2(3,3) = m1(1,1)*m1(2,2)-m1(1,2)*m1(2,1)

detm = m1(1,1)*m1(2,2)*m1(3,3) + m1(1,2)*m1(2,3)*m1(3,1) &
     + m1(1,3)*m1(2,1)*m1(3,2) - m1(1,3)*m1(2,2)*m1(3,1) &
     - m1(1,2)*m1(2,1)*m1(3,3) - m1(1,1)*m1(2,3)*m1(3,2) 

m2(:,:) = m2(:,:)/detm

end subroutine

!--------------------------------------------------------------------------------------------------------------
function l2g(atype)
implicit none
!convert Local ID to Global ID 
!--------------------------------------------------------------------------------------------------------------
real(8),intent(IN) :: atype
integer :: l2g,ity

ity = nint(atype)
l2g = nint((atype-ity)*1d13)

return
end function

!--------------------------------------------------------------------------------------------------------------
subroutine xu2xs(mcx, rreal, rnorm, nmax)
! update normalized coordinate from real coordinate. Subtract obox to make them local. 
use md_context
!--------------------------------------------------------------------------------------------------------------
implicit none

type(md_context_type),intent(in) :: mcx
real(8),intent(in),allocatable :: rreal(:,:)
real(8),intent(inout),allocatable :: rnorm(:,:)
integer,intent(in) :: nmax

integer :: i
real(8) :: rr(3)

do i=1,nmax
   rr(1:3) = rreal(i,1:3)
   rnorm(i,1)=sum(mcx%HHi(1,1:3)*rr(1:3))
   rnorm(i,2)=sum(mcx%HHi(2,1:3)*rr(1:3))
   rnorm(i,3)=sum(mcx%HHi(3,1:3)*rr(1:3))
   rnorm(i,1:3) = rnorm(i,1:3) - mcx%OBOX(1:3)
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine xs2xu(mcx,rnorm,rreal,nmax)
! update real coordinate from normalized coordinate
use md_context
!--------------------------------------------------------------------------------------------------------------

type(md_context_type),intent(in) :: mcx
real(8),intent(in),allocatable :: rnorm(:,:)
real(8),intent(inout),allocatable :: rreal(:,:)
integer,intent(in) :: nmax

real(8) :: rr(3)

do i=1,nmax 
   rr(1:3) = rnorm(i,1:3) + mcx%OBOX(1:3)
   rreal(i,1)=sum(mcx%HH(1,1:3,0)*rr(1:3))
   rreal(i,2)=sum(mcx%HH(2,1:3,0)*rr(1:3))
   rreal(i,3)=sum(mcx%HH(3,1:3,0)*rr(1:3))
enddo

end subroutine

!-----------------------------------------------------------------------
subroutine ScaleTemperature(mcx, ffp, mpt, treq, atype, v)
use md_context; use ff_params; use mpi_vars
!-----------------------------------------------------------------------
implicit none

type(md_context_type),intent(in) :: mcx
type(forcefield_params),intent(in) :: ffp
type(mpi_var_type),intent(in) :: mpt
real(8),intent(in) :: treq

real(8),allocatable :: atype(:), v(:,:)

integer :: i,ity
real(8) :: Ekinetic, ctmp

do i=1,mcx%NATOMS
   ity=nint(atype(i))
   Ekinetic=0.5d0*ffp%mass(ity)*sum(v(i,1:3)*v(i,1:3))
   ctmp = (treq*UTEMP0)/( Ekinetic*UTEMP )
   v(i,1:3)=sqrt(ctmp)*v(i,1:3)
enddo

call LinearMomentum(mcx, ffp, mpt, atype, v)

return
end

!-----------------------------------------------------------------------
subroutine LinearMomentum(mcx, ffp, mpt, atype, v)
use md_context; use ff_params; use mpi_vars
!-----------------------------------------------------------------------
implicit none

type(md_context_type),intent(in) :: mcx
type(mpi_var_type),intent(in) :: mpt
type(forcefield_params),intent(in) :: ffp
real(8),allocatable :: atype(:), v(:,:)

integer :: i,ity,ierr
real(8) :: mm,vCM(3),sbuf(4),rbuf(4)

!--- get the local momentum and mass.
vCM(:)=0.d0;  mm = 0.d0
do i=1,mcx%NATOMS
   ity = nint(atype(i))
   vCM(1:3)=vCM(1:3) + ffp%mass(ity)*v(i,1:3)
   mm = mm + ffp%mass(ity)
enddo

sbuf(1)=mm; sbuf(2:4)=vCM(1:3)
call MPI_ALLREDUCE(sbuf, rbuf, size(sbuf), MPI_DOUBLE_PRECISION, MPI_SUM, mpt%mycomm, ierr)
mm=rbuf(1); vCM(1:3)=rbuf(2:4)

!--- get the global momentum
vCM(:)=vCM(:)/mm

!--- set the total momentum to be zero 
do i=1,mcx%NATOMS
   v(i,1:3) = v(i,1:3) - vCM(1:3)
enddo

return
end


!----------------------------------------------------------------
subroutine GetBoxParams(H,la,lb,lc,angle1,angle2,angle3)
!----------------------------------------------------------------
implicit none
real(8),intent(inout) :: H(3,3)
real(8),intent(in) :: la,lb,lc, angle1,angle2,angle3
real(8) :: hh1, hh2 , lal, lbe, lga
real(8) :: pi=atan(1.d0)*4.d0

!--- convet unit for angles
lal=angle1*pi/180.d0
lbe=angle2*pi/180.d0
lga=angle3*pi/180.d0

!--- construct H-matrix
hh1=lc*(cos(lal)-cos(lbe)*cos(lga))/sin(lga)
hh2=lc*sqrt( 1.d0-cos(lal)**2-cos(lbe)**2-cos(lga)**2 + &
             2*cos(lal)*cos(lbe)*cos(lga) )/sin(lga)

H(1,1)=la;          H(2,1)=0.d0;        H(3,1)=0.d0
H(1,2)=lb*cos(lga); H(2,2)=lb*sin(lga); H(3,2)=0.d0
H(1,3)=lc*cos(lbe); H(2,3)=hh1;         H(3,3)=hh2

return
end subroutine

!----------------------------------------------------------------
subroutine UpdateBoxParams(mcx, rxp)
use md_context; use rxmd_params; use ff_params
!----------------------------------------------------------------
implicit none

type(md_context_type),intent(inout) :: mcx
type(rxmd_param_type),intent(in) :: rxp

!--- get volume 
mcx%MDBOX = &
mcx%HH(1,1,0)*(mcx%HH(2,2,0)*mcx%HH(3,3,0) - mcx%HH(3,2,0)*mcx%HH(2,3,0)) + &
mcx%HH(2,1,0)*(mcx%HH(3,2,0)*mcx%HH(1,3,0) - mcx%HH(1,2,0)*mcx%HH(3,3,0)) + &
mcx%HH(3,1,0)*(mcx%HH(1,2,0)*mcx%HH(2,3,0) - mcx%HH(2,2,0)*mcx%HH(1,3,0))

!--- get inverse of H-matrix
call matinv(mcx%HH,mcx%HHi)

!--- local box dimensions (a temporary use of lbox)
mcx%LBOX(1)=mcx%lata/rxp%vprocs(1)
mcx%LBOX(2)=mcx%latb/rxp%vprocs(2)
mcx%LBOX(3)=mcx%latc/rxp%vprocs(3)

!--- get the number of linkedlist cell per domain
mcx%cc(1:3)=int(mcx%LBOX(1:3)/mcx%maxrc)

!--- local system size in the unscaled coordinate.
mcx%LBOX(1:3) = 1.d0/rxp%vprocs(1:3)

!--- get the linkedlist cell dimensions (normalized)
mcx%lcsize(1:3) = mcx%LBOX(1:3)/mcx%cc(1:3)

!--- get origin of local MD box in the scaled coordiate.
mcx%OBOX(1:3) = mcx%LBOX(1:3)*mcx%vID(1:3)

return
end

end module
