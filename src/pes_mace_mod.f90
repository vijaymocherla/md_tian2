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

    type emt_mace
      private
      character(len = 256) :: model_path
      !type(MaceModel) :: calc
    end type emt_mace

    type(emt_mace) :: pes_mace

    ! pbc conditions are hard coded here
    contains 

    subroutine read_mace(atoms, inp_unit)

      use run_config, only : simparams

      type(universe), intent(inout) :: atoms
      integer, intent(in) :: inp_unit

      integer :: nwords, ios = 0
      character(len=max_string_length) :: buffer
      character(len=max_string_length) :: words(100)
      integer  :: idx1, idx2, ntypes, param_counter
      character(len=*), parameter :: err = "Error in read_mace: "

      ntypes = simparams%nprojectiles+simparams%nlattices
      pes_mace%model_path = default_string

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
      type(MaceModel) :: model
      ! allocate node energy
      allocate(node_energy(atoms%natoms))
      ! pass the model path
      print *, "Calculating MACE forces"
      ! initialize the model
      pes_mace%calc = MaceModel(adjustl(trim(pes_mace%model_path)))      
      call pes_mace%calc%print()
      call pes_mace%calc%calculate(.true., atoms%natoms, atoms%simbox, pbc, atoms%atomic_numbers, atoms%r, total_energy, node_energy, atoms%f, virial)
      atoms%epot = total_energy
      call pes_mace%calc%deallocate()
      deallocate(node_energy)
    end subroutine compute_mace


end module pes_mace_mod
