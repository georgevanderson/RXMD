!-------------------------------------------------------------------------------------------
module ff_params
!-------------------------------------------------------------------------------------------
!  This module stores the parameters input from the rxmda.in file.  This file is 
!    intentionally very similar to the input setup used by Adri van Duin at Caltech (to aid
!    in cross-checking and updates).  However, please note that there are a few changes made
!    aside from the obvious f77 -> f90 switch.  I have tried to clearly note these. 
!-------------------------------------------------------------------------------------------
!Taper function 
real(8),parameter :: rctap0 = 10.d0 ![A]

! hydrogen bonding interaction cutoff
real(8),parameter :: rchb = 10.d0 ![A]
real(8),parameter :: rchb2 = rchb*rchb

!Coulomb Energy (eq. 22)  
real(8),parameter:: Cclmb0 = 332.0638d0 ! [kcal/mol/A] line 2481 in poten.f
real(8),parameter:: Cclmb0_qeq = 14.4d0 ! [ev]
real(8),parameter:: CEchrge = 23.02d0   ! [ev]
real(8) :: Cclmb = Cclmb0

!Independant Parameters

type forcefield_params 

   integer :: nso    !Number of different types of atoms
   integer :: nboty  !Number of different bonds given
   
   ! Atom Dependant (ie where they appear in input file - not implementation in code)
   character(2),allocatable :: atmname(:)      !holds the Chemical Abbrev for each atomtype
   real(8),allocatable :: Val(:),Valboc(:)  !Valency of atomtype (norm, boc) 
   real(8),allocatable :: mass(:)           !mass of atomtype
   
   real(8),allocatable :: pbo1(:), pbo2(:), pbo3(:)   !Bond Order terms
   real(8),allocatable :: pbo4(:), pbo5(:), pbo6(:)   !Bond Order terms
   real(8),allocatable :: pboc1(:), pboc2(:), pboc3(:)  !Bond Order correction terms (f1-5)
   real(8),allocatable :: pboc4(:), pboc5(:)            !Bond Order correction terms (f1-5)  
   real(8),allocatable :: v13cor(:) !<kn>
   
   real(8),allocatable :: rat(:),rapt(:),vnq(:)   !r0s/r0p/r0pp for like bonds 
   real(8),allocatable :: ovc(:)  !a flag to apply fn4 and fn5 !<kn>
   
   real(8),allocatable :: Desig(:),Depi(:),Depipi(:)  !Bond Energy parameters (eq. 6)
   real(8),allocatable :: pbe1(:),pbe2(:)             !Bond Energy parameters (eq. 6) 
   
   real(8),allocatable :: Vale(:)                      !Lone Pair Energy parameters (eq. 7)
   real(8),allocatable :: plp1(:), nlpopt(:), plp2(:)  !Lone Pair Energy parameters (eq.8-10)  
   
   real(8),allocatable :: povun1(:), povun2(:), povun3(:), povun4(:)   !Overcoordination Energy (eq. 11)
   real(8),allocatable :: povun5(:), povun6(:), povun7(:), povun8(:)   !Undercoordination Energy (eq. 12)
   
   real(8),allocatable :: pval1(:), pval2(:), pval3(:), pval4(:), pval5(:)   !Valency Angle Energy (eq. 13a-g)
   real(8),allocatable :: pval6(:), pval7(:), pval8(:), pval9(:), pval10(:)   
   real(8),allocatable :: Valangle(:), theta00(:)
   
   real(8),allocatable :: ppen1(:), ppen2(:), ppen3(:), ppen4(:)   !Penalty Energy (eq. 14ab)
   
   real(8),allocatable :: pcoa1(:), pcoa2(:), pcoa3(:), pcoa4(:)   !Conjugation (3 body) Energy (eq.15)
   real(8),allocatable :: Valval(:)
   
   real(8),allocatable :: ptor1(:), ptor2(:), ptor3(:), ptor4(:)   !Torsional Energy Terms (eq.16abc)
   real(8),allocatable :: V1(:), V2(:), V3(:)
   
   real(8),allocatable :: pcot1(:), pcot2(:)  !Conjugation (4body) Energy (eq. 17ab)
   
   real(8),allocatable :: phb1(:), phb2(:), phb3(:), r0hb(:)   !Hydrogren Bond Energy (eq. 18)
   
   real(8),allocatable :: Dij(:,:), alpij(:,:), rvdW(:,:), gamW(:,:)  !Van der Waals Energy (eq. 21ab)
   real(8) :: pvdW1, pvdW1h, pvdW1inv
   
   !Taper function 
   real(8) :: rctap, rctap2, CTap(0:7)
   
   real(8),allocatable :: gam(:), gamij(:,:)  
   
   !Charge Equilibration part, <chi>  electronegativity   <eta> stiffness
   real(8),allocatable :: chi(:), eta(:)
   
   !End Lost Parameters Listing
   
   !Not Understood Parameters: 
   real(8),allocatable :: bom(:)
   
   ! 2-atom combo dependant: 
   real(8),allocatable :: r0s(:,:)       !Bond Order terms
   real(8),allocatable :: r0p(:,:)       !  "" 
   real(8),allocatable :: r0pp(:,:)      !Bond Order terms 
   
   ! <inxn2> is type of 2-body interaction 1=C-C, 2=H-C, 3=H-H
   integer,allocatable :: inxn2(:,:)
   
   integer,allocatable :: inxn3(:,:,:)
   integer,allocatable :: inxn3hb(:,:,:)
   integer,allocatable :: inxn4(:,:,:,:)
   
   !Saved calculations to prevent having to recalc lots of times.
   real(8),allocatable :: cBOp1(:), cBOp3(:), cBOp5(:)
   real(8),allocatable :: pbo2h(:), pbo4h(:), pbo6h(:)  
   
   !NOTE: for debugging purpose variables
   real(8)  :: vpar30,vpar1,vpar2
   
   !--- <switch> flag to omit pi and double pi bond in bond-order prime calculation.
   real(8),allocatable :: switch(:,:) 

end type forcefield_params

contains 

!-------------------------------------------------------------------------------------------
SUBROUTINE GETPARAMS(fp, ffFileName, ffFileHeader)
!-------------------------------------------------------------------------------------------
!  This subroutine is designed solely to obtain the parameters used in the Ecalc.f90 
!    program from the rxmd.in input file.  It is similar to the input routine used
!    by Adri in "reac.f::ffinput"
!-------------------------------------------------------------------------------------------
implicit none

type(forcefield_params) :: fp

real(8),parameter :: pi=3.14159265358979d0

character(*),intent(in) :: ffFileName  ! force field parm file
character(*),intent(inout) :: ffFileHeader ! 1st line of the FF file

integer :: i,j,inxn   !counters for initialization
integer :: i0,i1,i2,i3,i4,ih  !Counters: # corresp to #-atom depend 

!--- Parameters that count number of entries in each field:
!integer :: nso     !* of atom types (in mod) 
!integer :: nboty  !# of bond terms given (in mod)
integer :: nodmty  !# of off-diag terms given   
integer :: npar, nvaty, ntoty, nhbty

!--- Readin Fields converted to Parameters:
integer :: nodm1, nodm2  !ID bonds to alter in off-diag terms
real(8) :: deodmh,rodmh,godmh ! off-diag term of Evdw parameters
real(8) :: rsig,rpi,rpi2 !temp storage for r0s, etc terms in off-diag part
integer :: typea,typeb   !Temp storage for filling inxn2(:,:) table

!--- NULL Transfer Fields not needed in program (used to calc other values): 
real(8) :: dnull
real(8),allocatable :: vpar(:), bo131(:), bo132(:), bo133(:)
real(8),allocatable :: rvdw1(:), eps(:), alf(:), vop(:)

dnull = 0.d0
!--- Start Getting Parameters
open(4,file=trim(adjustl(ffFileName)),status="old")
read(4,'(a100)') ffFileHeader

read(4,*) npar  !num of parameters (independ of atom choice)

allocate(vpar(npar))

do i0=1, npar 
   read(4,1300) vpar(i0)  !temp variable...some apparently depend on atype
enddo
  
!--- Constant parameters reset to actual use:
fp%pvdW1 = vpar(29)
fp%pvdW1h = 0.5d0*fp%pvdW1 
fp%pvdW1inv = 1.d0/fp%pvdW1

!--- a small modification in sigma-bond prime <kn>
fp%vpar30 = vpar(30)  

!--- Parameters with 1-atom Depend,  Nr of types of atoms
read(4,'(i3)') fp%nso    

!--- Allocation of variables:
allocate(fp%rat(fp%nso),fp%rapt(fp%nso),fp%vnq(fp%nso))
allocate(fp%r0s(fp%nso,fp%nso),fp%r0p(fp%nso,fp%nso),fp%r0pp(fp%nso,fp%nso))
allocate(fp%Val(fp%nso),fp%Valboc(fp%nso))
allocate(fp%mass(fp%nso))
allocate(bo131(fp%nso),bo132(fp%nso),bo133(fp%nso))
allocate(fp%inxn2(fp%nso,fp%nso),fp%inxn3(fp%nso,fp%nso,fp%nso))
allocate(fp%inxn3hb(fp%nso,fp%nso,fp%nso),fp%inxn4(fp%nso,fp%nso,fp%nso,fp%nso))
allocate(fp%atmname(fp%nso))
allocate(fp%Vale(fp%nso), fp%plp1(fp%nso), fp%nlpopt(fp%nso), fp%plp2(fp%nso))
allocate(fp%povun2(fp%nso),fp%povun3(fp%nso),fp%povun4(fp%nso))
allocate(fp%povun5(fp%nso),fp%povun6(fp%nso),fp%povun7(fp%nso),fp%povun8(fp%nso))
!--- Valency Terms (j-dependancy only):
allocate(fp%pval3(fp%nso),fp%pval5(fp%nso), fp%Valangle(fp%nso),fp%Valval(fp%nso))
!--- Van der Waals Terms:
allocate(rvdw1(fp%nso), fp%rvdW(fp%nso,fp%nso), eps(fp%nso), fp%Dij(fp%nso,fp%nso))
allocate(alf(fp%nso), fp%alpij(fp%nso,fp%nso))
allocate(vop(fp%nso), fp%gamW(fp%nso,fp%nso))

!--- Coulomb & Charge equilibration:
allocate(fp%chi(fp%nso), fp%eta(fp%nso), fp%gam(fp%nso), fp%gamij(fp%nso,fp%nso))
 
!--- Parameters that still don't depend on atom type yet
fp%plp1(1:fp%nso) = vpar(16)
fp%povun3(1:fp%nso) = vpar(33)
fp%povun4(1:fp%nso) = vpar(32)
fp%povun6(1:fp%nso) = vpar(7)
fp%povun7(1:fp%nso) = vpar(9)
fp%povun8(1:fp%nso) = vpar(10)

!--- skip 3 lines
read(4,*)
read(4,*)
read(4,*)

do i1=1, fp%nso  !collect info on each type of atom
   read(4,1200) fp%atmname(i1),fp%rat(i1),fp%Val(i1),fp%mass(i1),rvdw1(i1),eps(i1),fp%gam(i1),fp%rapt(i1),fp%Vale(i1)
   read(4,1250) alf(i1),vop(i1),fp%Valboc(i1),fp%povun5(i1),dnull,fp%chi(i1),fp%eta(i1),dnull
   read(4,1250) fp%vnq(i1),fp%plp2(i1),dnull,bo131(i1),bo132(i1),bo133(i1),dnull,dnull   
   read(4,1250) fp%povun2(i1), fp%pval3(i1),dnull,fp%Valval(i1),fp%pval5(i1)
enddo

fp%nlpopt(1:fp%nso) = 0.5d0*(fp%Vale(1:fp%nso) - fp%Val(1:fp%nso))
!--- duplicate values
fp%Valangle(1:fp%nso) = fp%Valboc(1:fp%nso)

!--- Calc default r0s, r0p, r0pp:
do i1=1,fp%nso
   do i2=1,fp%nso
!--- Terms for the Bond Order Calculation:
      fp%r0s(i1,i2) = 0.5d0*(fp%rat(i1)+fp%rat(i2))
      fp%r0p(i1,i2) = 0.5d0*(fp%rapt(i1)+fp%rapt(i2))
      fp%r0pp(i1,i2) = 0.5d0*(fp%vnq(i1)+fp%vnq(i2))   
    
!--- Terms used in van der Waals calc: 
      fp%rvdW(i1,i2) = sqrt( 4.d0*rvdw1(i1)*rvdw1(i2) )
      fp%Dij(i1,i2) = sqrt( eps(i1)*eps(i2) )
      fp%alpij(i1,i2) = sqrt( alf(i1)*alf(i2) )
      fp%gamW(i1,i2) = sqrt( vop(i1)*vop(i2) )  
      fp%gamij(i1,i2) = ( fp%gam(i1)*fp%gam(i2) )**(-1.5d0) !<- gamcco in reac.f
   enddo
enddo  

!--- Parameters with 2-atom Depend:
read(4,1100) fp%nboty  !# of bonds' params given 

!--- Allocation of variables:
allocate(fp%pbo1(fp%nboty),fp%pbo2(fp%nboty),fp%pbo3(fp%nboty))
allocate(fp%pbo4(fp%nboty),fp%pbo5(fp%nboty),fp%pbo6(fp%nboty),fp%bom(fp%nboty))
allocate(fp%pboc1(fp%nboty),fp%pboc2(fp%nboty),fp%pboc3(fp%nboty),fp%pboc4(fp%nboty),fp%pboc5(fp%nboty))
allocate(fp%desig(fp%nboty),fp%depi(fp%nboty),fp%depipi(fp%nboty),fp%pbe1(fp%nboty),fp%pbe2(fp%nboty)) 
allocate(fp%povun1(fp%nboty),fp%ovc(fp%nboty),fp%v13cor(fp%nboty))

!--- skip one line
read(4,*)

fp%inxn2(1:fp%nso,1:fp%nso) = 0   !allows later flag to tell when a combination is not allowed
ih=0
do i2=1,fp%nboty
    ih=ih+1
    read(4,1400) typea,typeb,fp%Desig(ih),fp%Depi(ih),fp%Depipi(ih),fp%pbe1(ih),fp%pbo5(ih),fp%v13cor(ih),fp%pbo6(ih),fp%povun1(ih)
    read(4,1450) fp%pbe2(ih),fp%pbo3(ih),fp%pbo4(ih),fp%bom(ih),fp%pbo1(ih),fp%pbo2(ih),fp%ovc(ih),dnull  
    fp%inxn2(typea,typeb) = ih
    fp%inxn2(typeb,typea) = ih
enddo


!--- TEMP (required by input file backwards setup)
fp%pboc1(1:fp%nboty) = vpar(1)
fp%pboc2(1:fp%nboty) = vpar(2)

!--- for debugging <kn>
fp%vpar1 = vpar(1)
fp%vpar2 = vpar(2)
 
do i1=1,fp%nso
  do i2=1,fp%nso 
    inxn = fp%inxn2(i1,i2)
    if(inxn/=0) then 
       fp%pboc3(inxn) = sqrt(bo132(i1)*bo132(i2)) ! be careful about variable name
       fp%pboc4(inxn) = sqrt(bo131(i1)*bo131(i2)) ! bo132 -> pboc2, bo131->pbo4  
       fp%pboc5(inxn) = sqrt(bo133(i1)*bo133(i2))
      endif
  enddo
enddo


!--- Changes to off-diagonal terms:
read(4,1100) nodmty  !# of off-diag terms 
do i2=1, nodmty
   read(4,1400) nodm1,nodm2,deodmh,rodmh,godmh,rsig,rpi,rpi2
   if(rsig.GT.0.d0) fp%r0s(nodm1,nodm2)=rsig
   if(rsig.GT.0.d0) fp%r0s(nodm2,nodm1)=rsig 
   if(rpi.GT.0.d0)  fp%r0p(nodm1,nodm2)=rpi
   if(rpi.GT.0.d0)  fp%r0p(nodm2,nodm1)=rpi
   if(rpi2.GT.0.d0) fp%r0pp(nodm1,nodm2)=rpi2
   if(rpi2.GT.0.d0) fp%r0pp(nodm2,nodm1)=rpi2
   if (rodmh.GT.0.d0) fp%rvdW(nodm1,nodm2)=2.0*rodmh 
   if (rodmh.GT.0.d0) fp%rvdW(nodm2,nodm1)=2.0*rodmh
   if (deodmh.GT.0.d0) fp%Dij(nodm1,nodm2)=deodmh
   if (deodmh.GT.0.d0) fp%Dij(nodm2,nodm1)=deodmh
   if (godmh.GT.0.d0) fp%alpij(nodm1,nodm2)=godmh
   if (godmh.GT.0.d0) fp%alpij(nodm2,nodm1)=godmh
enddo

!!--- Derived 2body parameters 
allocate(fp%cBOp1(fp%nboty), fp%cBOp3(fp%nboty), fp%cBOp5(fp%nboty))
allocate(fp%pbo2h(fp%nboty), fp%pbo4h(fp%nboty), fp%pbo6h(fp%nboty))

!--- <switch> flag to omit pi and double pi bond.
allocate(fp%switch(1:3,fp%nboty))

fp%switch(:,:)=0
do i=1,fp%nso
   do j=1,fp%nso
   inxn = fp%inxn2(i,j)

   if(inxn/=0) then

!!--- In BOp calculation, <switch> will be multiplied to <BOp> to remove
!!--- BOpi and BOpipi for bonding interaction of atoms with a hydrogen.
       if((fp%rat(i)>0.d0)  .and. fp%rat(j)>0.d0 )  fp%switch(1,inxn)=1
       if((fp%rapt(i)>0.d0) .and. fp%rapt(j)>0.d0 ) fp%switch(2,inxn)=1
       if((fp%vnq(i)>0.d0)  .and. fp%vnq(j)>0.d0 )  fp%switch(3,inxn)=1

      if(fp%r0s(i,j)<=0.d0) then
         fp%cBOp1(inxn) = 0.d0
      else
         fp%cBOp1(inxn) = fp%pbo1(inxn)/(fp%r0s(i,j)**fp%pbo2(inxn))
      endif
      if(fp%r0p(i,j)<=0.d0) then
         fp%cBOp3(inxn) = 0.d0
      else
         fp%cBOp3(inxn) = fp%pbo3(inxn)/(fp%r0p(i,j)**fp%pbo4(inxn))
      endif

      if(fp%r0pp(i,j)<=0.d0) then
         fp%cBOp5(inxn) = 0.d0
      else
         fp%cBOp5(inxn) = fp%pbo5(inxn)/(fp%r0pp(i,j)**fp%pbo6(inxn))
      endif

      fp%pbo2h(inxn) = 0.5d0*fp%pbo2(inxn)
      fp%pbo4h(inxn) = 0.5d0*fp%pbo4(inxn)
      fp%pbo6h(inxn) = 0.5d0*fp%pbo6(inxn)
   endif
   enddo
enddo

!--- Input Valency Terms from Input File
fp%inxn3(1:fp%nso,1:fp%nso,1:fp%nso) = 0
read(4,1100) nvaty

allocate(fp%pval1(nvaty),fp%pval2(nvaty),fp%pval4(nvaty),fp%pval6(nvaty))
allocate(fp%pval7(nvaty),fp%pval8(nvaty),fp%pval9(nvaty),fp%pval10(nvaty))
allocate(fp%theta00(nvaty))
allocate(fp%ppen1(nvaty),fp%ppen2(nvaty),fp%ppen3(nvaty),fp%ppen4(nvaty))
allocate(fp%pcoa1(nvaty),fp%pcoa2(nvaty),fp%pcoa3(nvaty),fp%pcoa4(nvaty))

do i=1, nvaty
   read(4,1500) i1,i2,i3,fp%theta00(i),fp%pval1(i),fp%pval2(i),fp%pcoa1(i),fp%pval7(i),fp%ppen1(i),fp%pval4(i)
   fp%inxn3(i1,i2,i3) = i
   fp%inxn3(i3,i2,i1) = i
enddo

fp%inxn3(1,2,2)=0 !react.f, line 4933
fp%inxn3(2,2,1)=0 !react.f, line 4933

!--- Valency Terms which do not depend on inxn type:
fp%pval6(1:nvaty) = vpar(15)
fp%pval8(1:nvaty) = vpar(34)
fp%pval9(1:nvaty) = vpar(17)
fp%pval10(1:nvaty) = vpar(18)
!--- Penalty Terms which do not depend on inxn type:
fp%ppen2(1:nvaty) = vpar(20)
fp%ppen3(1:nvaty) = vpar(21)
fp%ppen4(1:nvaty) = vpar(22)
!--- 3body Conjugation Terms which do not depend on type:
fp%pcoa2(1:nvaty) = vpar(3)
fp%pcoa3(1:nvaty) = vpar(39)
fp%pcoa4(1:nvaty) = vpar(31)
!--- theta00 given in degrees, but used in radians. Convert by:
fp%theta00(1:nvaty) = (pi/180.d0)*fp%theta00(1:nvaty) 


read(4,1100) ntoty
allocate(fp%ptor1(ntoty),fp%ptor2(ntoty),fp%ptor3(ntoty),fp%ptor4(ntoty))
allocate(fp%V1(ntoty),fp%V2(ntoty),fp%V3(ntoty))
allocate(fp%pcot1(ntoty),fp%pcot2(ntoty))

fp%inxn4(1:fp%nso,1:fp%nso,1:fp%nso,1:fp%nso) = 0  
do i=1,ntoty
   read(4,1600)i1,i2,i3,i4,fp%V1(i),fp%V2(i),fp%V3(i),fp%ptor1(i),fp%pcot1(i),dnull,dnull 
!--- Set up inxn4 lookup reference array
      if(i1==0) then   !condensed input, means that all i1,i4 for this arrangement of i2,i3 are the same
         do i1=1,fp%nso
         do i4=1,fp%nso
         if(fp%inxn4(i1,i2,i3,i4)==0.and.fp%inxn4(i1,i3,i2,i4)==0) then
            fp%inxn4(i1,i2,i3,i4) = i
            fp%inxn4(i4,i2,i3,i1) = i
            fp%inxn4(i1,i3,i2,i4) = i
            fp%inxn4(i4,i3,i2,i1) = i
          endif
          enddo
          enddo
      else
          fp%inxn4(i1,i2,i3,i4) = i
          fp%inxn4(i4,i2,i3,i1) = i
          fp%inxn4(i1,i3,i2,i4) = i
          fp%inxn4(i4,i3,i2,i1) = i 
      endif
enddo

!and a few which don't depend on type 
fp%ptor2(1:ntoty) = vpar(24)
fp%ptor3(1:ntoty) = vpar(25)
fp%ptor4(1:ntoty) = vpar(26)  
fp%pcot2(1:ntoty) = vpar(28)

!--- Input Hydrogen Bond Terms
fp%inxn3hb(1:fp%nso,1:fp%nso,1:fp%nso) = 0
read(4,1100) nhbty
allocate(fp%phb1(nhbty),fp%phb2(nhbty),fp%phb3(nhbty),fp%r0hb(nhbty))

do i=1,nhbty
   read(4,1500) i1,i2,i3,fp%r0hb(i),fp%phb1(i),fp%phb2(i),fp%phb3(i)
   fp%inxn3hb(i1,i2,i3) = i    !Note: inxn3hb(i,j,k) /= inxn3hb(k,j,i)
enddo


!--- close parameter file "ffield"
close(4)

!--- Formats:
1100 format (i3,2x,a2,3x,3d22.15)
1200 format (1x,a2,10f9.4)
1250 format (3x,10f9.4)
1300 format (f10.4)
1400 format (2i3,8f9.4)
1450 format (6x,8f9.4)
1500 format (3i3,7f9.4)
1600 format (4i3,7f9.4)


!--- coefficient of coulomb energy
!Cclmb = Cclmb0  !Eclmb

!--- Parameters for charge variable routine. 
!--- In original parameter, chiEEM and etaEEM are given in [ev], not [kcal/mol]
!--- Definition of the stiffness parameter <eta> is different from 
!--- the original code and our code. It's need to be multiplied by 2.
fp%eta(:) = fp%eta(:)*2.d0

fp%rctap = rctap0
fp%rctap2 = fp%rctap**2

fp%CTap(0:7)=(/1.d0, 0.d0, 0.d0, 0.d0,   -35.d0/(fp%rctap)**4, &
          84.d0/(fp%rctap)**5, -70.d0/(fp%rctap)**6, &
          20.d0/(fp%rctap)**7 /)

END SUBROUTINE

end module
