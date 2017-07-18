!-------------------------------------------------------------------------------------------
module atom_vars
!-------------------------------------------------------------------------------------------
! position, atom type, velocity, force & charge
type atom_var_type
  real(8),allocatable,dimension(:) :: atype, q
  real(8),allocatable,dimension(:,:) :: pos, v, f
end type atom_var_type

end module
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
module mpi_vars
!-------------------------------------------------------------------------------------------
include 'mpif.h'

type mpi_var_type
   integer :: myid, nprocs, mycomm
   integer :: provided
   integer :: ierr
end type mpi_var_type

contains

subroutine GetMPIVariables(mpt)
   implicit none

   type(mpi_var_type) ::mpt

   mpt%mycomm = MPI_COMM_WORLD

   call MPI_INIT(mpt%ierr)
   !call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,mpt%provided,ierr)
   call MPI_COMM_RANK(mpt%mycomm, mpt%myid, mpt%ierr)
   call MPI_COMM_SIZE(mpt%mycomm, mpt%nprocs, mpt%ierr)

end subroutine GetMPIVariables

end module mpi_vars

!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
module cmdline_args
!-------------------------------------------------------------------------------------------

integer,parameter :: MAXPATHLENGTH=256

type cmdline_arg_type

  logical :: saveRunProfile=.false.
  logical :: isFF=.false., isData=.false., isMDparm=.false.
  character(MAXPATHLENGTH) :: FFPath="ffield", DataDir="DAT", ParmPath="rxmd.in"

end type cmdline_arg_type 

contains 

subroutine GetCmdLineArgs(cla)
implicit none

type(cmdline_arg_type) :: cla
integer :: i
character(MAXPATHLENGTH) :: argv

!--- read FF file, output dir, MD parameter file paths from command line
do i=1, command_argument_count()
   call get_command_argument(i,argv)
   select case(adjustl(argv))
     case("--help","-h")
       print'(a)', "--ffield ffield --outDir DAT --rxmdin rxmd.in"
       stop
     case("--ffield", "-ff")
       call get_command_argument(i+1,argv)
       cla%FFPath=adjustl(argv)
     case("--outDir", "-o")
       call get_command_argument(i+1,argv)
       cla%DataDir=adjustl(argv)
     case("--rxmdin", "-in")
       call get_command_argument(i+1,argv)
       cla%ParmPath=adjustl(argv)
     case("--profile")
       cla%saveRunProfile=.true.
     case default
   end select

enddo
end subroutine

end module cmdline_args

!-------------------------------------------------------------------------------------------
module rxmd_params
!-------------------------------------------------------------------------------------------

type rxmd_param_type
! MD mode
   integer :: mdmode

! one time step in femto second.
   real(8) :: dt
   
! number of MD steps for current run. 
   integer :: ntime_step 
   
! requested temperature. 
   real(8) :: treq 
   
! velocity scaling factor.
   real(8) :: vsfact
   
! temperature scaling steps.
   integer :: sstep

! call OUTPUT() every fstep MD steps and save data into files if their flags are on.
   integer :: fstep

! call PRINTE() every pstep MD steps to show value of energy terms.
   integer :: pstep
   
! number of processers in [xyz] directions. 
   integer :: vprocs(3)
   
!--- QEq related variables. 
! flag to run QEq routine: 0-No QEq, 1-CG, 2-Extended Lagrangian
   integer :: isQEq
   
! Number of MAXimum iteration in QEq routine. 
   integer :: NMAXQEq
   
! energy convergence criterion in QEq routine.
   real(8) :: QEq_tol

!  call QEq() every qstep MD steps.
   integer :: qstep

!  variables for extended Lagurangian method.
   real(8) :: Lex_fqs=1.0, Lex_k=2.d0

! flags to tell which file to be saved. 
   logical :: isBinary, isBondFile, isPDB

! tolerance of energy convergence in conjugate gradient.
   real(8) :: ftol   

end type rxmd_param_type

contains

!-------------------------------------------------------------------------------------------
subroutine GetRxmdParams(rxp, paramPath)
implicit none
!-------------------------------------------------------------------------------------------

type(rxmd_param_type),intent(out) :: rxp
character(*),intent(in) :: paramPath

!--- read MD control parameters
open(1, file=trim(paramPath), status="old")
read(1,*) rxp%mdmode
read(1,*) rxp%dt, rxp%ntime_step
read(1,*) rxp%treq, rxp%vsfact, rxp%sstep
read(1,*) rxp%fstep, rxp%pstep
read(1,*) rxp%vprocs(1:3)
read(1,*) rxp%isQEq, rxp%NMAXQEq, rxp%QEq_tol, rxp%qstep
read(1,*) rxp%Lex_fqs, rxp%Lex_k
read(1,*) rxp%isBinary, rxp%isBondFile, rxp%isPDB
read(1,*) rxp%ftol
close(1)

end subroutine 

end module rxmd_params


!-------------------------------------------------------------------------------------------
module atoms
!-------------------------------------------------------------------------------------------

integer,parameter :: MAXPATHLENGTH=256

character(MAXPATHLENGTH) :: RunProfilePath="profile.dat"
integer,parameter :: RunProfileFD=30 ! file descriptor for summary file

!--- For array size statistics
!  1-NATOMS, 2-nbrlist, 3-nbrlist for qeq, 4-NBUFFER for move, 5-NBUFFER for copy
!  6-NBUFFER for qeq
integer,parameter :: nmaxas=5
integer,allocatable :: maxas(:,:)

!--- lattice parameters 
real(8) :: lata,latb,latc,lalpha,lbeta,lgamma

integer :: ierr, myparity(3), vID(3)

!<NE_COPY>,<NE_MOVE>,<NE_CPBK> :: Number of Elements to COPY, MOVE atoms and CoPy BacK force. 
integer,parameter :: MODE_COPY = 1, MODE_MOVE = 2, MODE_CPBK = 3
integer,parameter :: MODE_QCOPY1 = 4, MODE_QCOPY2 = 5, MODE_STRESSCALC = 6

integer,parameter :: NE_COPY = 10, NE_MOVE = 12
integer,parameter :: NE_QCOPY1 = 2, NE_QCOPY2 = 3, NE_STRESSCALC = 6

#ifdef STRESS
integer,parameter :: NE_CPBK = 10
#else
integer,parameter :: NE_CPBK = 4
#endif

!<MAXLAYERS> MAXimum # of linkedlist cell LAYERS.
integer,parameter :: MAXLAYERS=5
integer,parameter :: MAXLAYERS_NB=10
        
! <target_node> stores partner node ID in the 6-communications. 
! if targe_node(i)==-1, the node doesn't have a partner in i-direction.
integer :: target_node(6)

real(8),allocatable:: rc(:), rc2(:)   !<RCUT>: cutoff length for sigma-bonding.
real(8),allocatable:: rcpi(:), rcpp(:)!      : cutoff length for other bonding.

real(8),parameter :: MINBOSIG = 1d-3      !<minBOsig>: criterion to decide <rc> 
real(8),parameter :: MINBO0 = 1d-4       !<minBO0>: cutoff bond order 
real(8),parameter :: cutof2_esub = 1d-4
!real(8),parameter :: cutof2_bo = 1.d-2
real(8),parameter :: cutof2_bo = 1.d-3
integer,parameter :: is_idEh = 1

!real(8),parameter :: MINBOSIG = 1d-4      !<minBOsig>: criterion to decide <rc> 
!real(8),parameter :: MINBO0 = 0.d0       !<minBO0>: cutoff bond order 
!real(8),parameter :: cutof2_esub = 0.d0
!real(8),parameter :: cutof2_bo = 1d-4
!integer,parameter :: is_idEh = 0

!integer :: NBUFFER=5000
!integer,parameter :: MAXNEIGHBS=50  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
!integer,parameter :: MAXNEIGHBS10=200 !<MAXNEIGHBS>: Max # of Ngbs within 10[A]. 

integer :: NBUFFER=30000
integer,parameter :: MAXNEIGHBS=30  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
integer,parameter :: MAXNEIGHBS10=700 !<MAXNEIGHBS>: Max # of Ngbs within 10[A]. 

integer,parameter :: NMINCELL=3  !<NMINCELL>: Nr of minimum linkedlist cell <-> minimum grain size.
real(8),parameter :: MAXANGLE= 0.999999999999d0 
real(8),parameter :: MINANGLE=-0.999999999999d0
real(8),parameter :: NSMALL = 1.d-10
real(8) :: maxrc                        !<maxRCUT>: Max cutoff length. used to decide lcsize.

real(8),parameter :: pi=3.14159265358979d0

! atomic stress tensor
real(8),allocatable :: astr(:,:) 
real(8) :: pint(3,3)

real(8) :: HH(3,3,0:1), HHi(3,3), MDBOX, LBOX(0:3), OBOX(1:3) !MD box, local MD box, origin of box.
integer :: NATOMS         !local # of atoms
integer(8) :: GNATOMS     !global # of atoms
integer :: ALLATOMS

!<llist> Linked List
!<header> header atom of linkedlist cell.
!<nacell> Nr of atoms in a likedlist cell.
integer,allocatable :: llist(:), header(:,:,:), nacell(:,:,:)

!<nbllist> Linked List for non-bonding interaction
!<nbheader> header atom of linkedlist cell for non-bonding interaction
!<nbnacell> Nr of atoms in a likedlist cell for non-bonding interaction
integer,allocatable :: nbllist(:), nbheader(:,:,:), nbnacell(:,:,:)

!<nbrlist> neighbor list, <nbrindx> neighbor index
integer,allocatable :: nbrlist(:,:), nbrindx(:,:)

!<nbplist> neighbor list of nonbonding interaction, non-bonding pair list
integer,allocatable :: nbplist(:,:)

!--- Passed between Elnpr and E3body
real(8),allocatable :: nlp(:), dDlp(:) !Number of Lone Pairs, its derivatives.
real(8),allocatable :: deltalp(:)

! TE: Total Energy,  KE: Kinetic Energy,  PE :: Potential Energies
!  0-Esystem, 1-Ebond, 2-Elp, 3-Eover, 4-Eunder, 5-Eval, 6-Epen
!  7-Ecoa,  8-Etors, 9-Econj, 10-Ehbond, 11-Evdwaals, 12-Ecoulomb 13-Echarge
real(8) :: TE, KE, PE(0:13)
real(8) :: GTE, GKE, GPE(0:13)

! Two vectors electrostatic energy minimization 
real(8),allocatable :: qs(:),qt(:),gs(:), gt(:), hs(:), ht(:), hshs(:), hsht(:)

!--- variables for extended Lagrangian method ---
!<Lex_fqs> fraction between two QEq vectors
!<Lex_w> spring constant
real(8),allocatable :: qsfp(:),qsfv(:),qtfp(:),qtfv(:) 
real(8),allocatable :: hessian(:,:)
real(8) :: Lex_w=1.d0, Lex_w2=1.d0
 
! <lcsize> Linked list Cell SIZE. <cc> Nr of likedlist cell in local node.
real(8) :: lcsize(3), nblcsize(3)
integer :: cc(3), nbcc(3)

integer :: nmesh, nbnmesh
integer,allocatable :: mesh(:,:), nbmesh(:,:)

!--- Unit convetors. In original ReaxFF program, the units of length is [A]
!--- energy is [kcal/mol] and mass is [amu] respectively.
!--- Most of numbers written here are obtained below site.
!--- http://physics.nist.gov/cuu/Constants/Table/allascii.txt
!--- Bohr length
real(8),parameter :: Lbohr_a  = 0.5291772108d0   ! [A]
real(8),parameter :: Lbohr_m  = 0.5291772108d-10 ! [m]

!--- Electron rest mass
real(8),parameter :: Merest_amu = 5.48580152d-4  ! [amu]
real(8),parameter :: Merest_kg  = 9.10938d-31    ! [kg]
!--- Energy in Hartree
real(8),parameter :: Ehrtre_km = 6.2751d2        ! [kcal/mol]
real(8),parameter :: Ehrtre_ev  = 27.2113845d0   ! [eV]
real(8),parameter :: Ehrtre_j = 4.35974417d-18   ! [J] 
real(8),parameter :: Eev_kcal = 23.060538d0      ! [kcal/mol]

real(8),parameter :: Ekcal_j = 6.95016611d-21  ! [J]

!--- Boltzmann Constant
real(8),parameter :: BLTZMN = 1.3806503d-23  ! [m^2 kg s^-2 K^-1 ] 

real(8),parameter :: UTEMP0 = 503.398008d0    ! Ekcal_j/BLZMN [K]
real(8),parameter :: UTEMP = UTEMP0*2.d0/3.d0 ! [K]
real(8),parameter :: USTRS = 6.94728103d0     ! [GPa]
real(8),parameter :: UDENS = 1.66053886d0     ! [g/cc]
real(8),parameter :: UTIME = 1.d3/20.455d0    ! 1 = 1/20.445[ps] = 48.88780[fs]

!<nstep_qeq> counter of iteration
integer :: nstep_qeq

!-- variables for timing
integer,parameter :: Ntimer=30
integer :: it_timer(Ntimer)=0, it_timer_max(Ntimer)=0, it_timer_min(Ntimer)=0

! <nstep> current MD step, 
! <current_step> will be used for subsequent runs.
integer :: nstep=0, current_step

!--- <frcindx> FoRCe INDeX. Index to return calculated force to original atoms.
integer,allocatable :: frcindx(:)
integer :: copyptr(0:6)

!--- stress components
real(8) :: xx,yy,zz,yz,zx,xy
!--- atomic stress index
integer :: ia,ja

!--- cutoff range calculation. 
integer(8),allocatable :: natoms_per_type(:)

!--- dthm=dt/(2*mass), hmas=mass/2
real(8),allocatable :: dthm(:), hmas(:)

!--- potential teble
integer,parameter :: NTABLE=5000
real(8),allocatable :: TBL_Eclmb(:,:,:), TBL_Evdw(:,:,:), TBL_Eclmb_QEq(:,:)
real(8) :: UDR, UDRi

integer(8),allocatable :: ibuf8(:)

!--- FF parameter description
character(MAXPATHLENGTH) :: FFDescript

contains

!-----------------------------------------------------------------------------------------------------------------------
character(len=256) function rankToString(irank)
!-----------------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(in) :: irank

write(rankToString,*) irank
rankToString = adjustl(rankToString)

end function rankToString

!-------------------------------------------------------------------------------------------
function GetFileNameBase(dataDir, nstep) result(fileNameBase)
!-------------------------------------------------------------------------------------------
integer,intent(in) :: nstep
character(MAXPATHLENGTH) :: dataDir, fileNameBase
character(9) :: a9

if(nstep>=0) then
  write(a9,'(i9.9)') nstep
  fileNameBase=trim(adjustl(dataDir))//"/"//a9
else
  fileNameBase=trim(adjustl(dataDir))//"/rxff"
endif

end function

end module atoms

!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
module MemoryAllocator
!-------------------------------------------------------------------------------------------
implicit none
integer :: totalMemory=0

contains 

subroutine AllocatorD1D(array, imin, imax)
  implicit none
  integer,intent(in) :: imin, imax
  real(8),allocatable,dimension(:) :: array
  integer :: status
  
  allocate(array(imin:imax), stat=status)
  totalMemory = totalMemory + size(array)*8

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorD1D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorD1D(array)
  implicit none
  real(8),allocatable,dimension(:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorD1D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorD2D(array, imin1, imax1, imin2, imax2) 
  implicit none
  integer,intent(in) :: imin1, imax1, imin2, imax2
  real(8),allocatable,dimension(:,:) :: array
  integer :: status
  
  allocate(array(imin1:imax1,imin2:imax2), stat=status)
  totalMemory = totalMemory + size(array)*8

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorD2D: totalMemory = ', totalMemory, status

  return 
end subroutine

subroutine DeallocatorD2D(array) 
  implicit none
  real(8),allocatable,dimension(:,:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorD2D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorD3D(array, imin1, imax1, imin2, imax2, imin3, imax3) 
  implicit none
  integer,intent(in) :: imin1, imax1, imin2, imax2, imin3, imax3
  real(8),allocatable,dimension(:,:,:) :: array
  integer :: status
  
  allocate(array(imin1:imax1,imin2:imax2,imin3:imax3), stat=status)
  totalMemory = totalMemory + size(array)*8

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorD3D: totalMemory = ', totalMemory, status

  return 
end subroutine

subroutine DeallocatorD3D(array) 
  implicit none
  real(8),allocatable,dimension(:,:,:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorD3D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI1D(array, imin, imax) 
  implicit none
  integer,intent(in) :: imin, imax
  integer,allocatable,dimension(:) :: array
  integer :: status
  
  allocate(array(imin:imax), stat=status)
  totalMemory = totalMemory + size(array)*4

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorI1D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI1D(array) 
  implicit none
  integer,allocatable,dimension(:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorI1D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI2D(array, imin1, imax1, imin2, imax2) 
  implicit none
  integer,intent(in) :: imin1, imax1, imin2, imax2
  integer,allocatable,dimension(:,:) :: array
  integer :: status
  
  allocate(array(imin1:imax1,imin2:imax2), stat=status)
  totalMemory = totalMemory + size(array)*4

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorI2D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI2D(array) 
  implicit none
  integer,allocatable,dimension(:,:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorI2D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI3D(array, imin1, imax1, imin2, imax2, imin3, imax3) 
  implicit none
  integer,intent(in) :: imin1, imax1, imin2, imax2, imin3, imax3
  integer,allocatable,dimension(:,:,:) :: array
  integer :: status
  
  allocate(array(imin1:imax1,imin2:imax2,imin3:imax3), stat=status)
  totalMemory = totalMemory + size(array)*4

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorI3D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI3D(array) 
  implicit none
  integer,allocatable,dimension(:,:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorI3D: totalMemory = ', totalMemory, status

  return
end subroutine 

integer function GetTotalMemory() 
  GetTotalMemory = totalMemory
  return
end function

end module MemoryAllocator
!-------------------------------------------------------------------------------------------
