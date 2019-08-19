program main
  ! Design of a WWW to only work on the square planar case
  ! THIS PROGRAM IS NOT KIND TO THE GENERAL USER
  ! There is heavy reliance on the global variables in the top section,
  ! thus most subroutines appear in this source file
  use maths, only : inv
  use types
  use random

  implicit none

  real(dp), allocatable :: coords(:,:)                      !List of coordinates
  real(dp), allocatable :: coords_curr(:,:)                 !List of retained coordinates
  integer, allocatable :: adj(:,:)                          !Adjacency matrix
  character(len=2), allocatable :: sym(:)                   !Atom symbols
  character(len=10), allocatable :: label(:)                !Atom labels
  real(dp), dimension(2,3) :: lattice                       !Lattice parameters
  character(len=80)::filename                               !Filename
  integer :: i,j,k                                          !Dummy ints
  integer,allocatable :: M_idx(:)                           !Indecies of all metal centers
  integer :: M1,M2                                          !Chosen metal index
  integer,dimension(4) :: b1_idx,b2_idx                     !List of bridging indecies
  integer :: bridge                                         !Chosen bridging index
  integer :: b1,b2                                          !Index of chosen bonds to swap
  integer :: N_accept                                       !Number of externally accepted moves
  integer :: N_accept_probe                                 !Number of accepted move in thermostat probe
  integer :: N_bonds
  logical :: is_fract                                       !Are coordinates fractional?
  real(dp), dimension(3,3) ::f2c, c2f                       !Fractional to cartesion conversion mat
  real(dp) :: v                                             !Parameter for conversion mats
  real(dp) :: halfL(3)                                      !Half box lengths
  real(dp) :: Ef                                            !Final energies for internal
  real(dp) :: E_curr                                        !Current energy of box
  real(dp) :: timing(2)                                     !CPU timing
  real(dp) :: T                                             !Temperature
  real(dp) :: dT                                            !Temperature change
  real(dp) :: r                                             !Random number
  real(dp) :: inAccept                                      !interal acceptance

  ! Numbers to included in an input file--------------------------------------------------!
  real(dp) :: max_dist                                      !max dist for bonding
  integer :: N_x,N_y,N_z                                    !Supercell dimensions
  integer :: n_atoms                                        !number of atoms in unit cell
  integer :: in_moves                                       !Number of SWEEPS for internal MC
  real(dp) :: in_T0                                         !Internal starting temperature
  real(dp) :: in_Tf                                         !Internal final temp
  integer :: out_moves                                      !Moves for WWW MC
  real(dp) :: out_T0                                        !WWW starting temperature
  real(dp) :: out_Tf                                        !WWW final temperature
  real(dp) :: max_move                                      !Maximum move in angstroms
  real(dp) :: k_12                                          !Bond strength
  real(dp) :: r_12                                          !Eq Bond length (expectation value of bl)
  real(dp) :: r_hs                                          !Cuttoff for VDW repulsion, squared later
  real(dp) :: k3_1                                          !Angle rigidity, type 1, metal centers
  real(dp) :: k3_2                                          !Angle rigidity, type 2, metal centers
  real(dp) :: k3_3                                          !Angle rigidity, type 3, bridge centers
  real(dp) :: k_hs                                          !Vanderwalls repulsion (previously hard spheres)
  real(dp) :: theta1                                        !Angle, type 1,180, changed to cos(theta)
  real(dp) :: theta2                                        !Angle, type 2,90,changed to cos(theta)
  real(dp) :: theta3                                        !Angle, type 3,changed to cos(theta)
  character(len=2):: sym1                                   !Symbols of relevance for metal center
  character(len=100):: struct_in                            !Input cif/pdb filename
  character(len=80) :: out_dir                              !Output directory
  logical :: rflex                                          !Using expectation value vs fixed value
  character(len=20) :: thermostat                           !Case selction for thermostat
  real(dp) :: T_decay                                       !Temperature decay (adaptive T)
  real(dp) :: T_boost                                       !Temperature boost (adaptive T)
  real(dp) :: acc_min                                       !Minimum acceptance (adaptive T)
  real(dp) :: acc_max                                       !Maximum acceptance (adaptive T)
  integer :: T_interval                                     !Number of steps between T check (adaptve T)
  ! -------------------------------------------------------------------------------------------!

  ! Initialization-----------------------------------------------------------------------------!
  call cpu_time(timing(1))
  i = command_argument_count()
  if (i<1) then
     write(*,*)'Please start the program by writing ./SquareWWW [input file name]'
     write(*,*)''
     stop
  else
     call get_command_argument(1,filename)
  end if
  call get_input(filename)
  call system('rm -rf '//trim(out_dir))
  call system('mkdir -p '//trim(out_dir))
  is_fract = .true.
  allocate(coords(3,N_x*N_y*N_z*n_atoms),adj(N_x*N_y*N_z*n_atoms,N_x*N_y*N_z*n_atoms), &
       sym(N_x*N_y*N_z*n_atoms),label(N_x*N_y*N_z*n_atoms))
  in_moves = in_moves * size(sym)
  call read_cif(struct_in)

  ! Creating conversion matricies
  f2c=0._dp
  c2f=0._dp
  lattice(2,:) = lattice(2,:) * PI/180._dp
  v= sqrt(1-cos(lattice(2,1))**2-cos(lattice(2,2))**2 - cos(lattice(2,3))**2 + &
       2*cos(lattice(2,1))*cos(lattice(2,2))*cos(lattice(2,3)))
  f2c(1,1) = lattice(1,1)
  f2c(1,2) = lattice(1,2)*cos(lattice(2,3))
  f2c(2,2) = lattice(1,2)*sin(lattice(2,3))
  f2c(1,3) = lattice(1,3)*cos(lattice(2,2))
  f2c(2,3) = lattice(1,3)*(cos(lattice(2,1))-cos(lattice(2,2))*cos(lattice(2,3)))/sin(lattice(2,3))
  f2c(3,3) = lattice(1,3)*v/sin(lattice(2,3))
  c2f = inv(f2c)
  lattice(2,:) = lattice(2,:) * 180._dp/PI
  ! -------------------------------------------------------------------------------------------!

  ! Derived parameters-------------------------------------------------------------------------!
  halfL(1) = 0.5_dp*lattice(1,1)
  halfL(2) = 0.5_dp*lattice(1,2)
  halfL(3) = 0.5_dp*lattice(1,3)
  theta1 = cos(pi/180.*theta1)
  theta2 = cos(pi/180.*theta2)
  theta3 = cos(pi/180.*theta3)
  r_hs = r_hs**2                !For ease of calc w/o square roots
  N_bonds = 0
  ! -------------------------------------------------------------------------------------------!

  call init_random_seed()
  call gen_adj()
  filename = trim(out_dir)//'start.pdb'
  call write_pdb(filename)
  filename = trim(out_dir)//'start_connective.pdb'
  call write_full_connect(filename)
  filename = trim(out_dir)//'start.cif'
  call write_cif(filename)

  ! Building set of metal indecies
  j=0
  do i=1,size(sym)
     if(sym(i)==sym1) j=j+1
  end do
  allocate(M_idx(j))
  j=0
  do i=1,size(sym)
     if(sym(i)==sym1) then
        j=j+1
        M_idx(j) = i
     end if
  end do

  ! WWW looping
  ! Choose metal,collect bonded,choose bridge,collect bonded,choose M2,collect bonded,choose b to swap
  ! Needs to be in cartesian coordinates
  if(is_fract) call fract2cart()
  T = out_T0
  dT = -1.*log(out_Tf/out_T0)/out_moves
  allocate(coords_curr(size(coords,1),size(coords,2)))
  coords_curr = coords
  N_accept = 0
  N_accept_probe=0
  E_curr = total_energy()

  write(filename,"('internal_',I3.3,'.csv')")0
  filename = trim(out_dir)//filename
  call internal_MC(Ef,inAccept,filename)
  write(filename,"('progress_',I3.3,'.pdb')")0
  filename = trim(out_dir)//filename
  call write_pdb(filename)
  write(filename,"('progress_connective_',I3.3,'.pdb')")0
  filename = trim(out_dir)//filename
  call write_full_connect(filename)
  E_curr = Ef
  coords_curr = coords

  filename = trim(out_dir)//'time_data.csv'
  call write_time(filename)

  WWW_Loop: do i=1,out_moves

     ! Pick metal center 1
     call random_number(r)
     k=ceiling(r*size(M_idx))
     M1=M_idx(k)

     ! Pick Bridge
     k=1
     do j=1,size(adj,1)
        if(adj(j,M1)>0)then
           b1_idx(k) = j
           k=k+1
        end if
     end do
     call random_number(r)
     k=ceiling(r*size(b1_idx))
     bridge = b1_idx(k)

     ! Find metal center 2
     k=1
     do j=1,size(adj,1)
        if(adj(j,bridge)>0)then
           b2_idx(k) = j
           k=k+1
        end if
     end do
     do
        call random_number(r)
        k=ceiling(r*2)
        if(b2_idx(k) /= M1) then
           M2 = b2_idx(k)
           exit
        else
           cycle
        end if
     end do

     ! Find bonding (already done for M1)
     k=1
     b2_idx=0
     do j=1,size(adj,1)
        if(adj(j,M2) >0) then
           b2_idx(k) = j
           k=k+1
        end if
     end do

     ! Choose bonds and make changes (symmetrically in adjacency matrix)
     do
        call random_number(r)
        k=ceiling(r*size(b1_idx))
        if(b1_idx(k) /=bridge .and. adj(M2,b1_idx(k))==0) then
           b1 = b1_idx(k)
           exit
        else
           cycle
        end if
     end do
     do
        call random_number(r)
        k=ceiling(r*size(b2_idx))
        if(b2_idx(k) /= bridge .and. sym(b2_idx(k))==sym(b1) .and. adj(M1,b2_idx(k))==0) then
           b2 = b2_idx(k)
           exit
        else
           cycle
        end if
     end do

     adj(M1,b1) = 0
     adj(b1,M1) = 0
     adj(M1,b2) = 1
     adj(b2,M1) = 1
     adj(M2,b2) = 0
     adj(b2,M2) = 0
     adj(M2,b1) = 1
     adj(b1,M2) = 1

     ! Check for rings
     if (ringSearch3()) then
        adj(M1,b1) = 1
        adj(b1,M1) = 1
        adj(M1,b2) = 0
        adj(b2,M1) = 0
        adj(M2,b2) = 1
        adj(b2,M2) = 1
        adj(M2,b1) = 0
        adj(b1,M2) = 0
        goto 15
     end if

     ! Acceptance criterion
     ! If accepted coords_curr = coords, else coords=coords_curr
     if(mod(i,out_moves/100) == 0) then
        write(filename,"('internal_',I3.3,'.csv')")i/(out_moves/100)
        filename = trim(out_dir)//filename
        call internal_MC(Ef,inAccept,filename)
     else
        call internal_MC(Ef,inAccept)
     end if
     call random_number(r)
     if (exp(-1.*(Ef-E_curr)/T)>=r) then
        N_accept = N_accept +1
        N_accept_probe= N_accept_probe +1
        coords_curr = coords
        E_curr = Ef
     else
        coords = coords_curr
        adj(M1,b1) = 1
        adj(b1,M1) = 1
        adj(M1,b2) = 0
        adj(b2,M1) = 0
        adj(M2,b2) = 1
        adj(b2,M2) = 1
        adj(M2,b1) = 0
        adj(b1,M2) = 0
     end if

15   select case(trim(thermostat))
     case('adaptive')
        if(mod(i,T_interval)==0) then
           if(N_accept_probe>=acc_max*T_interval) then
              T=T*T_decay
           else if(N_accept_probe<=acc_min*T_interval) then
              T = T*T_boost
           end if
           N_accept_probe=0
        end if
     case('annealing')
        T = T*(1-dT)
     case default
        T = T*(1-dT)
     end select
     
     if(mod(i,out_moves/100) == 0) then
        write(filename,"('progress_',I3.3,'.pdb')")i/(out_moves/100)
        filename = trim(out_dir)//filename
        call write_pdb(filename)
        write(filename,"('progress_connective_',I3.3,'.pdb')")i/(out_moves/100)
        filename = trim(out_dir)//filename
        call write_full_connect(filename)
        filename = trim(out_dir)//'time_data.csv'
        call write_time(filename,i)
     else if (mod(i,out_moves/5000) == 0) then
        filename = trim(out_dir)//'time_data.csv'
        call write_time(filename,i)
     end if

  end do WWW_Loop

  call cpu_time(timing(2))
  write(*,*)'Time',timing(2)-timing(1)
  filename = trim(out_dir)//'finish.pdb'
  call write_pdb(filename)
  filename = trim(out_dir)//'finish.cif'
  call write_cif(filename)
  filename = trim(out_dir)//'finalConnectivity.pdb'
  call write_full_connect(filename)

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
Contains

  subroutine internal_MC(Ef,acceptance,filename)
    ! THE EFFECTIVE ZERO TEMPERATURE IS 1e-7
    integer :: idx                                          !Index for change
    integer :: i                                            !Dummy
    real(dp) :: r                                           !Random number
    real(dp) :: move(3)                                     !random move
    real(dp) :: dE                                          !Energy change
    real(dp),intent(inout) :: Ef                            !Final energy
    real(dp) :: T                                           !Internal Temperature
    real(dp) :: dT                                          !Change in temperature
    real(dp) :: r_12_old                                    !Expectation value of bond length prior to change
    real(dp),intent(out),optional :: acceptance             !Acceptance
    character(len=*),intent(in),optional :: filename        !Filename for output
    integer :: n_accept                                     !Number of moves accepted
    integer :: unit
    N_accept=0
    Ef = total_energy()
    T = in_T0
    r_12_old = r_12
    dT = -1*log(in_Tf/T)/in_moves
    if (present(filename)) then
       open(newunit=unit, file=trim(filename), status="replace")
       write(unit,*)'Move,T,E,Acceptance,<r_12>'
       write(unit,"(I8,',',ES12.5,',',ES12.5,',',ES12.5,',',F8.3)")0,T,Ef,0.,r_12
       close(unit)
    end if
    if(rflex) then
       do i=1,in_moves
          call random_number(r)
          idx = ceiling(r*size(coords,2))
          call random_number(move)
          move = (move-0.5)*max_move
          dE = local_change_rflex(move,idx)
          call random_number(r)
          if (exp(-1.*dE/T)>=r) then
             N_accept = N_accept+1
             Ef = Ef + dE
             r_12_old = r_12
          else
             coords(:,idx) = coords(:,idx) - move
             r_12 = r_12_old
          end if
          if (mod(i,in_moves/1000) == 0 .and. present(filename)) then
             open(newunit=unit, file=trim(filename),status="old", position="append")
             write(unit,"(I8,',',ES12.5,',',ES12.5,',',ES12.5,',',F8.3)")i,T,Ef,real(N_accept)/i*100.,r_12
             close(unit)
          end if
          T = T*(1-dT)
       end do
    else
       do i=1,in_moves
          call random_number(r)
          idx = ceiling(r*size(coords,2))
          call random_number(move)
          move = (move-0.5)*max_move
          dE = local_change(move,idx)
          call random_number(r)
          if (exp(-1.*dE/T)>=r) then
             N_accept = N_accept+1
             Ef = Ef + dE
          else
             coords(:,idx) = coords(:,idx) - move
          end if
          if (mod(i,in_moves/1000) == 0 .and. present(filename)) then
             open(newunit=unit, file=trim(filename),status="old", position="append")
             write(unit,"(I8,',',ES12.5,',',ES12.5,',',ES12.5,',',F8.3)")i,T,Ef,real(N_accept)/i*100.,r_12
             close(unit)
          end if
          T = T*(1-dT)
       end do
    end if
    if (present(acceptance)) acceptance = real(N_accept)/in_moves*100
  end subroutine internal_MC
  ! -----------------------------------------------------------------------------

  function local_change(move,idx) result(r)
    ! Returns change in energy by single movement
    ! Actually enacts the movement
    ! Updates the expectation value of bond length, r_12
    real(dp), intent(in) :: move(3)                         !MC move
    integer, intent(in) :: idx                              !Atom index to move
    real(dp) :: v0                                          !Initial local potential
    real(dp) :: vf                                          !Final local potential
    real(dp) :: r                                           !Return change in energy
    integer :: i,j,k,l                                      !Dummy integers
    integer,allocatable :: con(:)                           !Connected indecies
    real(dp) ::dist(3)                                      !Distance for bonds
    real(dp),dimension(3) :: v1,v2                          !Vectors for angles
    integer,dimension(16) :: nearest2                       !Second order nearest neighbor(with redundancies, and 4 coordination)
    v0=0._dp
    vf=0._dp
    r=0._dp
    
    ! Building 2nd order neighbors list
    nearest2=0
    k=1
    i=idx
    do j=1,size(adj,1)
       if(adj(j,i)>0) then
          nearest2(k) = j
          k=k+1
       end if
    end do
    do l=1,count(nearest2/=0)
       do j=1,size(adj,1)
          if(adj(j,nearest2(l))>0) then
             nearest2(k)=j
             k=k+1
          end if
       end do
    end do
    ! Building indecies of connected atoms
    j=1
    do i=1,size(adj,1)
       if (adj(i,idx)>0) then
          j=j+1
       end if
    end do
    allocate(con(j))
    con=0
    con(1) = idx
    j=2
    do i=1,size(adj,1)
       if (adj(i,idx)>0) then
          con(j) = i
          j=j+1
       end if
    end do

    ! Bond length changes
    do i=1,size(coords,2)
       if (adj(idx,i) >0) then
          dist = abs(coords(:,idx) - coords(:,i))
          do k=1,3
             if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
          end do
          v0 = v0 + k_12*(r_12-norm2(dist))**2
       end if
    end do

    ! Bond angle changes
    ! First condition is for those centered at the central atoms
    ! Else condition is for those centered at the connecting atoms
    do i=1,size(con)
       if(sym(con(i)) == sym1) then
          do j=1,size(coords,2)-1
             if (adj(j,con(i))>0) then
                do k=j+1,size(coords,2)
                   if(adj(k,con(i))>0) then
                      v1 = coords(:,k) - coords(:,con(i))
                      v2 =  coords(:,j) - coords(:,con(i))
                      do l=1,3
                         if(v1(l)>halfL(l))then
                            v1(l) = halfL(l) - v1(l)
                         else if(v1(l)<-1*halfL(l)) then
                            v1(l) = -1*halfL(l) - v1(l)
                         end if
                      end do
                      do l=1,3
                         if(v2(l)>halfL(l))then
                            v2(l) = halfL(l) - v2(l)
                         else if(v2(l)<-1*halfL(l)) then
                            v2(l) = -1*halfL(l) - v2(l)
                         end if
                      end do
                      if(sym(j)==sym(k)) then
                         v0 = v0 + k3_1*(theta1 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                      else
                         v0 = v0 + k3_2*(theta2 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                      end if
                   end if
                end do
             end if
          end do
       else
          do j=1,size(coords,2)-1
             if (adj(j,con(i))>0) then
                do k=j+1,size(coords,2)
                   if(adj(k,con(i))>0) then
                      v1 = coords(:,k) - coords(:,con(i))
                      v2 =  coords(:,j) - coords(:,con(i))
                      do l=1,3
                         if(v1(l)>halfL(l))then
                            v1(l) = halfL(l) - v1(l)
                         else if(v1(l)<-1*halfL(l)) then
                            v1(l) = -1*halfL(l) - v1(l)
                         end if
                      end do
                      do l=1,3
                         if(v2(l)>halfL(l))then
                            v2(l) = halfL(l) - v2(l)
                         else if(v2(l)<-1*halfL(l)) then
                            v2(l) = -1*halfL(l) - v2(l)
                         end if
                      end do
                      v0 = v0 + k3_3*(theta3 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                   end if
                end do
             end if
          end do
       end if
    end do
    ! Hard spheres changes r^(-12) potential, but using r^2
    do i=1,size(coords,2)
       if(any(i==nearest2(:))) cycle
       dist = abs(coords(:,idx) - coords(:,i))
       do k=1,3
          if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
       end do
       if (dot_product(dist,dist) < r_hs) then
          v0 = v0 + k_hs*(dot_product(dist,dist)**(-6))
          !v0=v0+1.
       end if
    end do

    ! Acual Movement
    coords(:,idx) = coords(:,idx) + move
   
    ! Bond length changes
    do i=1,size(coords,2)
       if (adj(idx,i) >0) then
          dist = abs(coords(:,idx) - coords(:,i))
          do k=1,3
             if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
          end do
          vf = vf + k_12*(r_12-norm2(dist))**2
       end if
    end do
    ! Bond angle changes
    ! First condition is for those centered at the central atoms
    ! Else condition is for those centered at the connecting atoms
    do i=1,size(con)
       if(sym(con(i)) == sym1) then
          do j=1,size(coords,2)-1
             if (adj(j,con(i))>0) then
                do k=j+1,size(coords,2)
                   if(adj(k,con(i))>0) then
                      v1 = coords(:,k) - coords(:,con(i))
                      v2 =  coords(:,j) - coords(:,con(i))
                      do l=1,3
                         if(v1(l)>halfL(l))then
                            v1(l) = halfL(l) - v1(l)
                         else if(v1(l)<-1*halfL(l)) then
                            v1(l) = -1*halfL(l) - v1(l)
                         end if
                      end do
                      do l=1,3
                         if(v2(l)>halfL(l))then
                            v2(l) = halfL(l) - v2(l)
                         else if(v2(l)<-1*halfL(l)) then
                            v2(l) = -1*halfL(l) - v2(l)
                         end if
                      end do
                      if(sym(j)==sym(k)) then
                         vf = vf + k3_1*(theta1 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                      else
                         vf = vf + k3_2*(theta2 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                      end if
                   end if
                end do
             end if
          end do
       else
          do j=1,size(coords,2)-1
             if (adj(j,con(i))>0) then
                do k=j+1,size(coords,2)
                   if(adj(k,con(i))>0) then
                      v1 = coords(:,k) - coords(:,con(i))
                      v2 =  coords(:,j) - coords(:,con(i))
                      do l=1,3
                         if(v1(l)>halfL(l))then
                            v1(l) = halfL(l) - v1(l)
                         else if(v1(l)<-1*halfL(l)) then
                            v1(l) = -1*halfL(l) - v1(l)
                         end if
                      end do
                      do l=1,3
                         if(v2(l)>halfL(l))then
                            v2(l) = halfL(l) - v2(l)
                         else if(v2(l)<-1*halfL(l)) then
                            v2(l) = -1*halfL(l) - v2(l)
                         end if
                      end do
                      vf = vf + k3_3*(theta3 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                   end if
                end do
             end if
          end do
       end if
    end do
    ! Hard spheres changes r^(-12) potential, but using r^2
    do i=1,size(coords,2)
       if(any(i==nearest2(:))) cycle
       dist = abs(coords(:,idx) - coords(:,i))
       do k=1,3
          if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
       end do
       if (dot_product(dist,dist) < r_hs) then
          vf = vf + k_hs*(dot_product(dist,dist)**(-6))          
          !vf=vf+1.
       end if
    end do
    r=vf-v0
  end function local_change
  ! -----------------------------------------------------------------------------

  function total_energy() result(r)
    ! Function first updates the r_12 expectation value then calculates total energy
    ! Ths is called at the start of each MC relaxation, after each bond swap
    real(dp) :: r                                           !Result, energy
    real(dp) ::dist(3)                                      !Distance
    real(dp),dimension(3) :: v1,v2                          !Vectors for angles
    integer :: i,j,k,l
    integer,dimension(16) :: nearest2                       !Second order nearest neighbor(with redundancies, and 4 coordination)
    r=0._dp
    if (rflex) call r12_recalc()
    ! Running through bonds first with reference to central atom
    do i=1,size(coords,2)
       if(sym(i) /= sym1) cycle
       do j=1,size(coords,2)
          if(adj(j,i)>0) then
             dist = abs(coords(:,j) - coords(:,i))
             do k=1,3
                if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
             end do
             r = r+ k_12*(r_12-norm2(dist))**2
          end if
       end do
    end do
    ! Bond angles NOT centered at sym1 (O,S in NbO based network)
    do i=1,size(coords,2)
       if(sym(i) == sym1) cycle
       do j=1,size(coords,2)-1
          if(adj(j,i)>0) then
             do k=j+1, size(coords,2)
                if(adj(k,i)>0) then
                   v1 = coords(:,k) - coords(:,i)
                   v2 =  coords(:,j) - coords(:,i)
                   do l=1,3
                      if(v1(l)>halfL(l))then
                         v1(l) = halfL(l) - v1(l)
                      else if(v1(l)<-1*halfL(l)) then
                         v1(l) = -1*halfL(l) - v1(l)
                      end if
                   end do
                   do l=1,3
                      if(v2(l)>halfL(l))then
                         v2(l) = halfL(l) - v2(l)
                      else if(v2(l)<-1*halfL(l)) then
                         v2(l) = -1*halfL(l) - v2(l)
                      end if
                   end do
                   r = r + k3_3*(theta3 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                end if
             end do
          end if
       end do
    end do
    ! Bond angles centered at sym1 (Nb in NbO)
    do i=1,size(coords,2)
       if(sym(i) /= sym1) cycle
       do j=1,size(coords,2)-1
          if(adj(j,i)>0) then
             do k=j+1, size(coords,2)
                if(adj(k,i)>0) then
                   v1 = coords(:,k) - coords(:,i)
                   v2 =  coords(:,j) - coords(:,i)
                   do l=1,3
                      if(v1(l)>halfL(l))then
                         v1(l) = halfL(l) - v1(l)
                      else if(v1(l)<-1*halfL(l)) then
                         v1(l) = -1*halfL(l) - v1(l)
                      end if
                   end do
                   do l=1,3
                      if(v2(l)>halfL(l))then
                         v2(l) = halfL(l) - v2(l)
                      else if(v2(l)<-1*halfL(l)) then
                         v2(l) = -1*halfL(l) - v2(l)
                      end if
                   end do
                   if(sym(j)==sym(k)) then
                      r = r + k3_1*(theta1 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                   else
                      r = r + k3_2*(theta2 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                   end if
                end if
             end do
          end if
       end do
    end do
    ! Hard spheres energy r^(-12) potential, but using r^2
    do i=1,size(coords,2)
       ! Building 2nd order neighbors list
       nearest2=0
       k=1
       do j=1,size(adj,1)
          if(adj(j,i)>0) then
             nearest2(k) = j
             k=k+1
          end if
       end do
       do l=1,count(nearest2/=0)
          do j=1,size(adj,1)
             if(adj(j,nearest2(l))>0) then
                nearest2(k)=j
                k=k+1
             end if
          end do
       end do
       do j=i+1,size(coords,2)
          if(any(j==nearest2(:))) cycle
          dist = abs(coords(:,j) - coords(:,i))
          do l=1,3
             if(dist(l)>halfL(l)) dist(l) =  2*halfL(l)-dist(l)
          end do
          if (dot_product(dist,dist) < r_hs) then
             r = r + k_hs*(dot_product(dist,dist)**(-6))             
             !r=r+1
          end if
       end do
    end do
  end function total_energy
  ! -----------------------------------------------------------------------------

  function local_change_rflex(move,idx) result(r)
    ! Returns change in energy by single movement
    ! Actually enacts the movement
    ! Updates the expectation value of bond length, r_12
    real(dp), intent(in) :: move(3)                         !MC move
    integer, intent(in) :: idx                              !Atom index to move
    real(dp) :: v0                                          !Initial local potential
    real(dp) :: vf                                          !Final local potential
    real(dp) :: r                                           !Return change in energy
    integer :: i,j,k,l                                      !Dummy integers
    integer,allocatable :: con(:)                           !Connected indecies
    real(dp) ::dist(3)                                      !Distance for bonds
    real(dp),dimension(3) :: v1,v2                          !Vectors for angles
    integer,dimension(16) :: nearest2                       !Second order nearest neighbor(with redundancies, and 4 coordination)
    v0=0._dp
    vf=0._dp
    r=0._dp

    ! Building 2nd order neighbors list
    nearest2=0
    k=1
    i=idx
    do j=1,size(adj,1)
       if(adj(j,i)>0) then
          nearest2(k) = j
          k=k+1
       end if
    end do
    do l=1,count(nearest2/=0)
       do j=1,size(adj,1)
          if(adj(j,nearest2(l))>0) then
             nearest2(k)=j
             k=k+1
          end if
       end do
    end do
    ! Building indecies of connected atoms
    j=1
    do i=1,size(adj,1)
       if (adj(i,idx)>0) then
          j=j+1
       end if
    end do
    allocate(con(j))
    con=0
    con(1) = idx
    j=2
    do i=1,size(adj,1)
       if (adj(i,idx)>0) then
          con(j) = i
          j=j+1
       end if
    end do

    ! R_12 expectation value changes
    ! Adding all bond length changes to 'local' potential
    do i=1,size(coords,2)
       if(sym(i) /= sym1) cycle
       do j=1,size(coords,2)
          if(adj(j,i)>0) then
             dist = abs(coords(:,j) - coords(:,i))
             do k=1,3
                if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
             end do
             v0 = v0 + k_12*(r_12-norm2(dist))**2
          end if
       end do
    end do

    ! First half of r_12 update
    r_12 = r_12 * N_bonds
    do i=1,size(coords,2)
       if (adj(idx,i) >0) then
          dist = abs(coords(:,idx) - coords(:,i))
          do k=1,3
             if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
          end do
          r_12 = r_12 - norm2(dist)
       end if
    end do

    ! Bond angle changes
    ! First condition is for those centered at the central atoms
    ! Else condition is for those centered at the connecting atoms
    do i=1,size(con)
       if(sym(con(i)) == sym1) then
          do j=1,size(coords,2)-1
             if (adj(j,con(i))>0) then
                do k=j+1,size(coords,2)
                   if(adj(k,con(i))>0) then
                      v1 = coords(:,k) - coords(:,con(i))
                      v2 =  coords(:,j) - coords(:,con(i))
                      do l=1,3
                         if(v1(l)>halfL(l))then
                            v1(l) = halfL(l) - v1(l)
                         else if(v1(l)<-1*halfL(l)) then
                            v1(l) = -1*halfL(l) - v1(l)
                         end if
                      end do
                      do l=1,3
                         if(v2(l)>halfL(l))then
                            v2(l) = halfL(l) - v2(l)
                         else if(v2(l)<-1*halfL(l)) then
                            v2(l) = -1*halfL(l) - v2(l)
                         end if
                      end do
                      if(sym(j)==sym(k)) then
                         v0 = v0 + k3_1*(theta1 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                      else
                         v0 = v0 + k3_2*(theta2 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                      end if
                   end if
                end do
             end if
          end do
       else
          do j=1,size(coords,2)-1
             if (adj(j,con(i))>0) then
                do k=j+1,size(coords,2)
                   if(adj(k,con(i))>0) then
                      v1 = coords(:,k) - coords(:,con(i))
                      v2 =  coords(:,j) - coords(:,con(i))
                      do l=1,3
                         if(v1(l)>halfL(l))then
                            v1(l) = halfL(l) - v1(l)
                         else if(v1(l)<-1*halfL(l)) then
                            v1(l) = -1*halfL(l) - v1(l)
                         end if
                      end do
                      do l=1,3
                         if(v2(l)>halfL(l))then
                            v2(l) = halfL(l) - v2(l)
                         else if(v2(l)<-1*halfL(l)) then
                            v2(l) = -1*halfL(l) - v2(l)
                         end if
                      end do
                      v0 = v0 + k3_3*(theta3 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                   end if
                end do
             end if
          end do
       end if
    end do
    ! Hard spheres changes r^(-12) potential, but using r^2
    do i=1,size(coords,2)
       if(any(i==nearest2(:))) cycle
       dist = abs(coords(:,idx) - coords(:,i))
       do k=1,3
          if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
       end do
       if (dot_product(dist,dist) < r_hs) then
          v0 = v0 + k_hs*(dot_product(dist,dist)**(-6))
          !v0=v0+1.
       end if
    end do

    ! Acual Movement
    coords(:,idx) = coords(:,idx) + move


    ! R_12 expectation value changes
    ! Second half of r_12 update
    do i=1,size(coords,2)
       if (adj(idx,i) >0) then
          dist = abs(coords(:,idx) - coords(:,i))
          do k=1,3
             if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
          end do
          r_12 = r_12 + norm2(dist)
       end if
    end do
    r_12 = r_12 / N_bonds

    ! Updating 'local' potential with new r_12
    do i=1,size(coords,2)
       if(sym(i) /= sym1) cycle
       do j=1,size(coords,2)
          if(adj(j,i)>0) then
             dist = abs(coords(:,j) - coords(:,i))
             do k=1,3
                if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
             end do
             vf = vf + k_12*(r_12-norm2(dist))**2
          end if
       end do
    end do
    ! Bond angle changes
    ! First condition is for those centered at the central atoms
    ! Else condition is for those centered at the connecting atoms
    do i=1,size(con)
       if(sym(con(i)) == sym1) then
          do j=1,size(coords,2)-1
             if (adj(j,con(i))>0) then
                do k=j+1,size(coords,2)
                   if(adj(k,con(i))>0) then
                      v1 = coords(:,k) - coords(:,con(i))
                      v2 =  coords(:,j) - coords(:,con(i))
                      do l=1,3
                         if(v1(l)>halfL(l))then
                            v1(l) = halfL(l) - v1(l)
                         else if(v1(l)<-1*halfL(l)) then
                            v1(l) = -1*halfL(l) - v1(l)
                         end if
                      end do
                      do l=1,3
                         if(v2(l)>halfL(l))then
                            v2(l) = halfL(l) - v2(l)
                         else if(v2(l)<-1*halfL(l)) then
                            v2(l) = -1*halfL(l) - v2(l)
                         end if
                      end do
                      if(sym(j)==sym(k)) then
                         vf = vf + k3_1*(theta1 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                      else
                         vf = vf + k3_2*(theta2 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                      end if
                   end if
                end do
             end if
          end do
       else
          do j=1,size(coords,2)-1
             if (adj(j,con(i))>0) then
                do k=j+1,size(coords,2)
                   if(adj(k,con(i))>0) then
                      v1 = coords(:,k) - coords(:,con(i))
                      v2 =  coords(:,j) - coords(:,con(i))
                      do l=1,3
                         if(v1(l)>halfL(l))then
                            v1(l) = halfL(l) - v1(l)
                         else if(v1(l)<-1*halfL(l)) then
                            v1(l) = -1*halfL(l) - v1(l)
                         end if
                      end do
                      do l=1,3
                         if(v2(l)>halfL(l))then
                            v2(l) = halfL(l) - v2(l)
                         else if(v2(l)<-1*halfL(l)) then
                            v2(l) = -1*halfL(l) - v2(l)
                         end if
                      end do
                      vf = vf + k3_3*(theta3 - dot_product(v1,v2)/norm2(v1)/norm2(v2))**2
                   end if
                end do
             end if
          end do
       end if
    end do
    ! Hard spheres changes r^(-12) potential, but using r^2
    do i=1,size(coords,2)
       if(any(i==nearest2(:))) cycle
       dist = abs(coords(:,idx) - coords(:,i))
       do k=1,3
          if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
       end do
       if (dot_product(dist,dist) < r_hs) then
          vf = vf + k_hs*(dot_product(dist,dist)**(-6))          
          !vf=vf+1.
       end if
    end do
    r=vf-v0
  end function local_change_rflex
  ! -----------------------------------------------------------------------------

  subroutine r12_recalc()
    ! Calculates the expectation value of r_12 for use in the potential
    integer :: i,j
    real(dp) :: dist(3)
    if(N_bonds<=0) N_bonds = count_bonds()
    r_12 = 0
    ! Running through bonds first with reference to central atom
    do i=1,size(coords,2)
       if(sym(i) /= sym1) cycle
       do j=1,size(coords,2)
          if(adj(j,i)>0) then
             dist = abs(coords(:,j) - coords(:,i))
             do k=1,3
                if(dist(k)>halfL(k)) dist(k) =  2*halfL(k)-dist(k)
             end do
             r_12 = r_12 + norm2(dist)
          end if
       end do
    end do
    r_12 = r_12 / N_bonds
  end subroutine r12_recalc
  ! -----------------------------------------------------------------------------

  function count_bonds() result(r)
    ! Particularly useful for first run of r12_exp
    integer :: r
    integer :: i,j

    r=0
    do i=1,size(coords,2)
       if(sym(i) /= sym1) cycle
       do j=1,size(coords,2)
          if(adj(j,i)>0) then
             r=r+1
          end if
       end do
    end do
  end function count_bonds
  ! -----------------------------------------------------------------------------

  function ringSearch3() result(r)
    ! Searches for 3 membered rings of metals at M1 and M2
    ! This will fail at capturing 2 membered rings
    ! Returns true if rings exist
    logical :: r                                            !True for ring size 3
    integer :: Mn1(4)                                       !Metal neighbor indecies step 1
    integer :: Ol1(4)                                       !Oxygen linker
    integer :: Mn2(12)                                      !Metal neighbor indecies step 2
    integer :: i,j,k,l

    r=.false.
    Mn1=0
    Mn2=0
    Ol1=0
    k=1
    do i=1,size(adj,2)
       if(adj(i,M1)>0) then                                 !Finding linker
          do j=1,size(adj,2)
             if(adj(j,i)>0 .and. j/=M1) then                !Finding node on other side of linker
                Ol1(k) = i
                Mn1(k) = j
                k = k+1
             end if
          end do
       end if
    end do

    l=1
    do i=1,size(Mn1)
       do j=1,size(adj,2)
          if(adj(j,Mn1(i))>0 .and. j/=Ol1(i)) then          !Finding linkers for 1st neighborMn1(i)
             do k=1,size(adj,2)
                if(adj(k,j)>0 .and. k/=Mn1(i)) then         !Finding node other side
                   Mn2(l) = k
                   l=l+1
                end if
             end do
          end if
       end do
    end do

    l=1
    do i=1,size(Mn2)
       do j=1,size(adj,2)
          if(adj(j,Mn2(i))>0) then
             do k=1,size(adj,2)
                if(adj(k,j)>0 .and. k==M1) then
                   r = .true.
                   return
                end if
             end do
          end if
       end do
    end do

    ! Repeating for Metal 2
    Mn1=0
    Mn2=0
    Ol1=0
    k=1
    do i=1,size(adj,2)
       if(adj(i,M2)>0) then                                 !Finding linker
          do j=1,size(adj,2)
             if(adj(j,i)>0 .and. j/=M2) then                !Finding node on other side of linker
                Ol1(k) = i
                Mn1(k) = j
                k = k+1
             end if
          end do
       end if
    end do

    l=1
    do i=1,size(Mn1)
       do j=1,size(adj,2)
          if(adj(j,Mn1(i))>0 .and. j/=Ol1(i)) then          !Finding linkers for 1st neighborMn1(i)
             do k=1,size(adj,2)
                if(adj(k,j)>0 .and. k/=Mn1(i)) then         !Finding node other side
                   Mn2(l) = k
                   l=l+1
                end if
             end do
          end if
       end do
    end do

    l=1
    do i=1,size(Mn2)
       do j=1,size(adj,2)
          if(adj(j,Mn2(i))>0) then
             do k=1,size(adj,2)
                if(adj(k,j)>0 .and. k==M2) then
                   r = .true.
                   return
                end if
             end do
          end if
       end do
    end do
  end function ringSearch3
  ! -----------------------------------------------------------------------------
  subroutine gen_adj()
    integer :: i,j
    adj = 0
    do i=1,size(coords,2)
       do j=1,size(coords,2)
          if (i==j)cycle
          if(cart_distance(coords(:,j),coords(:,i))<max_dist) then
             adj(j,i) = adj(j,i)+1
          end if
       end do
    end do
  end subroutine gen_adj
  ! -----------------------------------------------------------------------------

  function cart_distance(c1,c2) result(r)
    real(dp), intent(in) :: c1(3),c2(3)                     !Coordinates of points
    real(dp) :: r
    real(dp) :: dist(3)
    integer :: i
    dist = c1(:)-c2(:)
    do i=1,3
       if(dist(i)>0.5) then
          dist(i)=dist(i)-1.0
       else if (dist(i)<-0.5) then
          dist(i)=dist(i)+1.0
       end if
    end do
    if (is_fract) then
       dist=matmul(f2c,dist)
    end if
    r=norm2(dist)
  end function cart_distance
  ! -----------------------------------------------------------------------------

  subroutine fract2cart()
    integer :: i
    if(.not. is_fract) return
    do i=1,size(coords,2)
       coords(:,i) = matmul(f2c,coords(:,i))
    end do
    is_fract = .false.
  end subroutine fract2cart
  subroutine cart2fract()
    integer :: i
    if(is_fract) return
    do i=1,size(coords,2)
       coords(:,i) = matmul(c2f,coords(:,i))
    end do
    is_fract = .true.
  end subroutine cart2fract
  ! -----------------------------------------------------------------------------

  subroutine read_cif(filename)
    ! Bastardized cif reader, will only work for P1 cifs in this format
    ! Makes suppercell lists from variables in main
    character(len=80) ::filename                            !Filename
    character(len=80) :: line
    character(len=80) :: target
    real :: dump
    integer :: unit
    integer :: info
    integer :: i,j

    open(newunit=unit, file=trim(filename), status='old',IOSTAT = info)
    if(info /= 0 ) then
       write(*,*)'Failure in bastardized cif reader. File not found'
       write(*,*) filename
    end if
    target = '_cell_length_a'
    do while(trim(line)/=trim(target))
       read(unit,*) line
    end do
    backspace(unit)
    do i=1,2
       do j=1,3
          read(unit,*) line, lattice(i,j)
       end do
    end do
    target = '_atom_site_fract_z'
    do while (trim(line)/=trim(target))
       read(unit,*) line
    end do

    do i=1,n_atoms
       read(unit,*) label(i),sym(i), dump, coords(:,i)
    end do
    close(unit)

    do i=1,N_x-1
       do j=1,n_atoms
          coords(:,(i)*n_atoms+j) = coords(:,j) + real((/i,0,0/))
          label((i)*n_atoms+j) = label(j)
          sym((i)*n_atoms+j) = sym(j)
       end do
    end do

    do i=1,N_y-1
       do j=1,n_atoms*N_x
          coords(:,(i)*n_atoms*N_x+j) = coords(:,j) + real((/0,i,0/))
          label((i)*n_atoms*N_x+j) = label(j)
          sym((i)*n_atoms*N_x+j) = sym(j)
       end do
    end do

    do i=1,N_z-1
       do j=1,n_atoms*N_x*N_y
          coords(:,i*n_atoms*N_x*N_y+j) = coords(:,j) + real((/0,0,i/))
          label(i*n_atoms*N_x*N_y+j) = label(j)
          sym(i*n_atoms*N_x*N_y+j) = sym(j)
       end do
    end do

    ! Change lattice parameters and divide by N
    lattice(1,1) = lattice(1,1)*N_x
    lattice(1,2) = lattice(1,2)*N_y
    lattice(1,3) = lattice(1,3)*N_z
    coords(1,:) = coords(1,:)/N_x
    coords(2,:) = coords(2,:)/N_y
    coords(3,:) = coords(3,:)/N_z
  end subroutine read_cif

  subroutine write_cif(filename)
    ! Writes generic P1 cif
    character(len=*), intent(in) :: filename
    integer :: i,unit

    if(.not. is_fract) call cart2fract()
    open(newunit=unit, file=trim(filename), status="replace")
    write(unit, "('data_cif')")
    write(unit,*)
    write(unit, "('_cell_length_a',T24,f9.5)") lattice(1,1)
    write(unit, "('_cell_length_b',T24,f9.5)") lattice(1,2)
    write(unit, "('_cell_length_c',T24,f9.5)") lattice(1,3)
    write(unit, "('_cell_angle_alpha',T24,f9.5)") lattice(2,1)
    write(unit, "('_cell_angle_beta',T24,f9.5)") lattice(2,2)
    write(unit, "('_cell_angle_gamma',T24,f9.5)") lattice(2,3)
    write(unit,*)
    write(unit, "('_symmetry_space_group_name_H-M', T48,'P1')")
    write(unit, "('_symmetry_Int_Tables_number',T48,'1')")
    write(unit, "('_symmetry_cell_setting',T48,'triclinic')")
    write(unit,*)
    write(unit,"('loop_')")
    write(unit,*)'_atom_site_label'
    write(unit,*)'_atom_site_type_symbol'
    write(unit,*)'_atom_site_fract_x'
    write(unit,*)'_atom_site_fract_y'
    write(unit,*)'_atom_site_fract_z'
    write(unit,*)'_atom_site_occupancy'

    do i=1,size(coords,2)
       write(unit,"(a10,T14,a2,T24,f8.5,T40,f8.5,T56,f8.5,T72,f8.4)")&
            label(i),sym(i), coords(:,i), 1.0000
    end do
    close(unit)
  end subroutine write_cif
  ! -----------------------------------------------------------------------------

  subroutine write_pdb(filename)
    character(len=*),intent(in) :: filename                 !Filename
    character(len=*),parameter :: fmt = "(A6,I5,A1,A2,A3,A3,A1,A1,I4,A1,A3,3F8.3,2F6.2,A10,A2,A2)"
    real(dp) :: dist(3)                                     !Distance for connectivity acrosss bounds
    integer :: unit
    integer :: i
    if(is_fract) call fract2cart()
    open(newunit=unit, file=trim(filename), status="replace")
    write(unit,"('HEADER    cheap fortran WWW generation')")
    do i=1,size(sym)
       write(unit,fmt)adjustl('ATOM  '),i,'',adjustr(sym(i)),'','MOL','','H',0,'','',coords(:,i),1.0,0.0,'',adjustr(sym(i)),'0 '
    end do
    do i=1,size(coords,2)
       write(unit,"('CONECT',I5)",ADVANCE='NO') i
       do j=1,size(coords,2)
          if (i==j) cycle
          dist=coords(:,i) - coords(:,j)
          dist = matmul(c2f,dist)
          if (adj(i,j)>0 .and. all(abs(dist)<=0.5_dp)) then
             write(unit,'(I5)', ADVANCE='NO') j
          end if
       end do
       write(unit,*)
    end do
    close(unit)
  end subroutine write_pdb
  ! -----------------------------------------------------------------------------

  subroutine write_full_connect(filename)
    character(len=*),intent(in) :: filename                 !Filename
    character(len=*),parameter :: fmt = "(A6,I5,A1,A2,A3,A3,A1,A1,I4,A1,A3,3F8.3,2F6.2,A10,A2,A2)"
    integer :: unit
    integer :: i
    if(is_fract) call fract2cart()
    open(newunit=unit, file=trim(filename), status="replace")
    write(unit,"('HEADER    cheap fortran WWW generation with full connectivity')")
    do i=1,size(sym)
       write(unit,fmt)adjustl('ATOM  '),i,'',adjustr(sym(i)),'','MOL','','H',0,'','',coords(:,i),1.0,0.0,'',adjustr(sym(i)),'0 '
    end do
    do i=1,size(coords,2)
       write(unit,"('CONECT',I5)",ADVANCE='NO') i
       do j=1,size(coords,2)
          if (i==j) cycle
          if (adj(i,j)>0) then
             write(unit,'(I5)', ADVANCE='NO') j
          end if
       end do
       write(unit,*)
    end do
    close(unit)
  end subroutine write_full_connect
  ! -----------------------------------------------------------------------------

  subroutine write_time(filename,move)
    character(len=*),intent(in) :: filename                 !Filename
    integer :: unit
    logical, save :: FirstCall = .true.
    integer,intent(in),optional :: move
    if(FirstCall) then
       open(newunit=unit, file=trim(filename), status="replace")
       write(unit,*)'Swap,CurrentE,CurrentT,OutAccept,InAccept%,<r_12>'
       write(unit,"(I8,',',ES12.5,',',ES12.5,',',I8,',',ES12.5,',',F8.3)") 0,E_curr,T,0,0.,r_12
       close(unit)
       FirstCall =.false.
    else
       open(newunit=unit, file=trim(filename), status="old",position="append")
       write(unit,"(I8,',',ES12.5,',',ES12.5,',',I8,',',ES12.5,',',F8.3)")&
            move,E_curr,T,N_accept, inAccept, r_12
       close(unit)
    end if
  end subroutine write_time
  ! -----------------------------------------------------------------------------

  subroutine get_input(filename)
    ! Input reader which sets global parameters marked up top
    character(len=*), intent(in) :: filename                !Filename
    integer :: unit                                         !File unit
    integer :: info                                         !Error information. 0 for success
    character(len=80) :: errmsg                             !Error message
    character(len=120) :: line                              !Read line
    integer :: i_line                                       !Line counter
    real :: power                                           !Power measure for log based inputs

    info = 0
    open(newunit=unit,file=trim(filename), status='old',IOSTAT = info,IOMSG=errmsg)

    ! Quitting immediately if input is unsuccessful
    if(info /=0) then
       write(*,*) errmsg
       stop
    end if

    N_x = 1
    N_y = 1
    N_z = 1
    n_atoms = 6
    max_dist = 2.5
    in_moves =int(5e5)
    in_T0 = (10._dp**(2))
    in_Tf = 10._dp**(-5)
    out_moves = int(1e5)
    out_T0 = 100.0
    out_Tf = 0.005
    max_move = 1.0
    r_12 = 2.105
    k_12 = 27.0
    k3_1 = 4.0
    k3_2 = 4.0
    k3_3 = 2.0
    k_hs = 100.
    r_hs = 2.0
    theta1 = 180._dp
    theta2 = 90._dp
    theta3 = 145._dp
    sym1 = 'Pd'
    rflex = .false.
    thermostat = 'annealing'
    T_decay = 0.75
    T_boost = 4.d0/3.d0
    acc_min = 0.01
    acc_max = 0.02
    T_interval = 500
    
    i_line=0
    do
       info=0
       read(unit,'(A120)',iostat = info) line
       i_line = i_line+1
       if (info /= 0) exit
       line = adjustl(line)                                 !Ignoring leading blanks
       if(line(1:1) == '!') cycle
       if (line(1:6) == 'struct') then
          read(line(7:120),'(A113)',iostat=info) struct_in
          if(info/=0) then
             write(*,*)'Bad filename for input structure'
             stop
          end if
          struct_in = trim(adjustl(struct_in))
          cycle
       else if (line(1:6) == 'outdir') then
          read(line(7:120),'(A113)',iostat=info) out_dir
          if(info/=0) then
             write(*,*)'Bad name for directory'
             stop
          end if
          out_dir = trim(adjustl(out_dir))
          cycle
       end if
       call lower_case(line)
       if(line(1:6) == 'natoms') then
          read(line(7:len_trim(line)),*,iostat=info,iomsg=errmsg) n_atoms
          if(info/=0) then
             write(*,*)'Bad value for number of atoms in the structure'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:2) == 'nx') then
          read(line(3:len_trim(line)),*,iostat=info,iomsg=errmsg) N_x
          if(info/=0) then
             write(*,*)'Bad value for number of cells in the supercell X direction'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:2) == 'ny') then
          read(line(3:len_trim(line)),*,iostat=info,iomsg=errmsg) N_y
          if(info/=0) then
             write(*,*)'Bad value for number of cells in the supercell Y direction'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:2) == 'nz') then
          read(line(3:len_trim(line)),*,iostat=info,iomsg=errmsg) N_z
          if(info/=0) then
             write(*,*)'Bad value for number of cells in the supercell Z direction'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:8) == 'max_dist') then
          read(line(9:len_trim(line)),*,iostat=info,iomsg=errmsg) max_dist
          if(info/=0) then
             write(*,*)'Bad value for maximum distance of bonding in the structure'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:8) == 'in_moves')then
          read(line(9:len_trim(line)),*,iostat=info,iomsg=errmsg) power
          if(info/=0) then
             write(*,*)'Bad value for internal moves'
             write(*,*)errmsg
             stop
          end if
          in_moves = int(10**power)
       else if(line(1:5) == 'in_t0')then
          read(line(6:len_trim(line)),*,iostat=info,iomsg=errmsg) power
          if(info/=0) then
             write(*,*)'Bad value for internal starting temperature'
             write(*,*)errmsg
             stop
          end if
          in_T0 = 10._dp**power
       else if(line(1:5) == 'in_tf')then
          read(line(6:len_trim(line)),*,iostat=info,iomsg=errmsg) power
          if(info/=0) then
             write(*,*)'Bad value for internal final temperature'
             write(*,*)errmsg
             stop
          end if
          in_Tf = 10._dp**power
       else if(line(1:8) == 'max_move')then
          read(line(9:len_trim(line)),*,iostat=info,iomsg=errmsg) max_move
          if(info/=0) then
             write(*,*)'Bad value for maximum move distacne'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:9) == 'out_moves')then
          read(line(10:len_trim(line)),*,iostat=info,iomsg=errmsg) power
          if(info/=0) then
             write(*,*)'Bad value for number of external WWW moves'
             write(*,*)errmsg
             stop
          end if
          out_moves = int(10**power)
       else if(line(1:6) == 'out_t0')then
          read(line(7:len_trim(line)),*,iostat=info,iomsg=errmsg) power
          if(info/=0) then
             write(*,*)'Bad value for outer starting temperature'
             write(*,*)errmsg
             stop
          end if
          out_T0 = 10._dp**power
       else if(line(1:6) == 'out_tf')then
          read(line(7:len_trim(line)),*,iostat=info,iomsg=errmsg) power
          if(info/=0) then
             write(*,*)'Bad value for outer final temperature'
             write(*,*)errmsg
             stop
          end if
          out_Tf = 10._dp**power
       else if(line(1:5) == 'rbond') then
          read(line(6:len_trim(line)),*,iostat=info,iomsg=errmsg) r_12
          if(info/=0) then
             write(*,*)'Bad value for bond radius'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:5) == 'kbond') then
          read(line(6:len_trim(line)),*,iostat=info,iomsg=errmsg) k_12
          if(info/=0) then
             write(*,*)'Bad value for bond strength'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:10) == 'kanglenode') then
          read(line(11:len_trim(line)),*,iostat=info,iomsg=errmsg) k3_1
          if(info/=0) then
             write(*,*)'Bad value for angle rigidity at node'
             write(*,*)errmsg
             stop
          end if
          k3_2=k3_1
       else if(line(1:10) == 'kanglelink') then
          read(line(11:len_trim(line)),*,iostat=info,iomsg=errmsg) k3_3
          if(info/=0) then
             write(*,*)'Bad value for angle rigidity at linker'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:8) == 'kspheres') then
          read(line(9:len_trim(line)),*,iostat=info,iomsg=errmsg) k_hs
          if(info/=0) then
             write(*,*)'Bad value for hard spheres potential'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:8) == 'rspheres') then
          read(line(9:len_trim(line)),*,iostat=info,iomsg=errmsg) r_hs
          if(info/=0) then
             write(*,*)'Bad value for hard spheres potential radius'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:10) == 'anglelink') then
          read(line(11:len_trim(line)),*,iostat=info,iomsg=errmsg) theta3
          if(info/=0) then
             write(*,*)'Bad value for angle  at linker'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:5) == 'rflex') then
          read(line(6:len_trim(line)),*,iostat=info,iomsg=errmsg) line
          if(info/=0) then
             write(*,*)'Bad input for rflex'
             write(*,*)errmsg
             stop
          end if
          if (trim(line)=='on') rflex = .true.
       else if(line(1:5) == 'node1') then
          read(line(6:len_trim(line)),*,iostat=info,iomsg=errmsg)sym1
          if(info/=0) then
             write(*,*)'Bad input for symbol of node 1'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:10) == 'thermostat') then
          read(line(11:len_trim(line)),*,iostat=info,iomsg=errmsg) thermostat
          if(info/=0) then
             write(*,*)'Bad input for thermostat'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:7) == 't_decay') then
          read(line(8:len_trim(line)),*,iostat=info,iomsg=errmsg) T_decay
          if(info/=0) then
             write(*,*)'Bad value for hard spheres potential'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:7) == 't_boost') then
          read(line(8:len_trim(line)),*,iostat=info,iomsg=errmsg) T_boost
          if(info/=0) then
             write(*,*)'Bad value for hard spheres potential'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:10) == 't_interval') then
          read(line(11:len_trim(line)),*,iostat=info,iomsg=errmsg) T_interval
          if(info/=0) then
             write(*,*)'Bad value for hard spheres potential'
             write(*,*)errmsg
             stop
          end if
       else if(line(1:13) == 'accept_bounds') then
          read(line(14:len_trim(line)),*,iostat=info,iomsg=errmsg) acc_min,acc_max
          if(info/=0) then
             write(*,*)'Bad value for hard spheres potential'
             write(*,*)errmsg
             stop
          end if
          acc_min=acc_min/100._dp
          acc_max=acc_max/100._dp
       else if(line(1:1) ==' ') then
          cycle
       else
          write(*,*)'Unrecognized input in input file line',i_line
       end if
    end do
  end subroutine get_input

  subroutine check_loc_con(move)
    ! Checks local metal connectivity to ensure 2 O and 2 S at each center
    integer :: i,j
    integer :: o,s
    integer, intent(in) :: move
    do i=1,size(adj,1)
       if(sym(i) /= sym1) cycle
       o=0;s=0;
       do j=1,size(adj,2)
          if(adj(j,i)>0) then
             if(trim(sym(j)) == 'O')then
                o = o+1
             else if(trim(sym(j)) == 'S') then
                s = s+1
             end if
          end if
       end do
       if(o/=2 .or. s /=2) then
          write(*,*)'Failure in keeping local connectivity at move',move
          write(*,*) i,sym(i),coords(:,i)
          stop
       end if
    end do
  end subroutine check_loc_con
end program main
