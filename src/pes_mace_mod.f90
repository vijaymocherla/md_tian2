module pes_mace_mod
    use universe_mod
    use useful_things
    use mace, only: MaceModel

    implicit none
      ! cell & atoms
    logical, dimension(3), parameter :: pbc = (/.true., .true., .true./)
    integer, allocatable :: atomic_numbers(:)

    real(dp), dimension(6) :: virial
      ! MACE model

    type mace_pes
      private
      character(len = 256) :: model_path
      integer, allocatable :: atomic_numbers(:)
      type(MaceModel) :: calc
    end type mace_pes

    type(mace_pes) :: pes_mace

    ! pbc conditions are hard coded here
    contains 

    subroutine read_mace(atoms, inp_unit)

      use run_config, only : simparams

      type(universe), intent(inout) :: atoms
      integer, intent(in) :: inp_unit

      integer :: nwords, ios = 0
      integer :: i, j, natoms
      character(len=max_string_length) :: buffer
      character(len=max_string_length) :: words(100)
      integer  :: idx1, idx2, ntypes, param_counter
      character(len=*), parameter :: err = "Error in read_mace: "
      ! set local variables
      ntypes = simparams%nprojectiles+simparams%nlattices
      pes_mace%model_path = default_string
      natoms = atoms%natoms

      if (.not. allocated(pes_mace%atomic_numbers)) then
        allocate(pes_mace%atomic_numbers(natoms))
      end if
      ! get atomic numbers
      do i=1, natoms
          do j=1, 118
              if (trim(atoms%name(i)) == trim(elements(j)%symbol)) then
                  pes_mace%atomic_numbers(i) = elements(j)%atomic_number
              end if
          end do
      end do
      ! line should read something like "H   H   proj    proj"
      read(inp_unit, '(A)', iostat=ios) buffer
      call split_string(buffer, words, nwords)

      if (nwords /= 4) stop err // "need four entries in interaction-defining lines"

      if (words(3) == "proj" .and. words(4) == "proj" .or. &
          words(3) == "proj" .and. words(4) == "latt" .or. &
          words(3) == "latt" .and. words(4) == "proj" .or. &
          words(3) == "latt" .and. words(4) == "latt") then

          idx1 = get_idx_from_name(atoms, words(1), is_proj=(words(3)=="proj"))
          idx2 = get_idx_from_name(atoms, words(2), is_proj=(words(4)=="proj"))

          if (atoms%pes(idx1,idx2) /= default_int) then
              print *, err // "pes already defined for atoms", words(1), words(3), words(2), words(4)
              stop
          end if

      else
          print *, err // "interaction must be defined via 'proj' and 'latt' keywords"
          stop
      end if

      ! set the pes type in the atoms object
      atoms%pes(idx1,idx2) = pes_id_mace
      atoms%pes(idx2,idx1) = pes_id_mace

      param_counter = 1
      do
          read(inp_unit, '(A)', iostat=ios) buffer
          call split_string(buffer, words, nwords)

          ! pes block terminated, exit
          if (nwords == 0 .or. ios /= 0) then
              exit

          ! something went wrong
          else if (nwords /= 2) then
              stop "Error in the PES file: PES parameters must consist of key value pairs. A parameter block must be terminated by a blank line."
          end if

          call lower_case(words(1))

          select case (words(1))

              case ('model_path')   ! A^-1
                  read(words(2), '(A)') pes_mace%model_path
              case default
                  print *, "Error in the PES file: unknown MACE MODEL", words(2)
                  stop

          end select
          param_counter = param_counter + 1
      end do

  end subroutine read_mace


    subroutine compute_mace(atoms, flag)
      type(universe), intent(inout) :: atoms
      integer, intent(in) :: flag
      real(dp), allocatable :: node_energy(:)
      real(dp):: virial(6)
      real(dp) :: total_energy
      type(MaceModel) :: calc
      ! allocate node energy
      allocate(node_energy(atoms%natoms))
      ! pass the model path
      print *, "Calculating MACE forces"
      ! initialize the model
      calc = MaceModel(adjustl(trim(pes_mace%model_path)))      
      call calc%print()
      call calc%calculate(.true., atoms%natoms, atoms%simbox, pbc, pes_mace%atomic_numbers, atoms%r, total_energy, node_energy, atoms%f, virial)
      atoms%epot = total_energy
      call calc%deallocate()
      deallocate(node_energy)
    end subroutine compute_mace


end module pes_mace_mod
