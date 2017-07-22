module qeq_vars

type qeq_var_type
! Two vectors electrostatic energy minimization 
  real(8),allocatable :: qs(:),qt(:),gs(:), gt(:), hs(:), ht(:), hshs(:), hsht(:)

!--- variables for extended Lagrangian method ---
!<Lex_fqs> fraction between two QEq vectors
!<Lex_w> spring constant
  real(8),allocatable :: qsfp(:),qsfv(:),qtfp(:),qtfv(:)
  real(8),allocatable :: hessian(:,:)
  real(8) :: Lex_w=1.d0, Lex_w2=1.d0

end type

end module
