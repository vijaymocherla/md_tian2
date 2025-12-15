!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2021 Daniel J. Auerbach, Svenja M. Janke, Marvin Kammler,
!               Sascha Kandratsenka, Sebastian Wille
! Dynamics at Surfaces Department
! MPI for Biophysical Chemistry Goettingen, Germany
! Georg-August-Universitaet Goettingen, Germany
!
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program. If not, see http://www.gnu.org/licenses.
!############################################################################

module constants

    use, intrinsic :: iso_fortran_env

    implicit none
    integer, parameter :: dp = REAL64
    
    type atomic_element
        character(len=3) :: symbol
        integer :: atomic_number
        real(dp) :: mass
    end type atomic_element
    

    ! Physical and mathematical constants
    real(dp), parameter :: sqrt2          = sqrt(2.0_dp)
    real(dp), parameter :: isqrt2         = 1.0_dp/sqrt2
    real(dp), parameter :: sqrt3          = sqrt(3.0_dp)
    real(dp), parameter :: pi             = acos(-1.0_dp)
    real(dp), parameter :: kB             = 8.61733238496e-5_dp       ! eV / K
    real(dp), parameter :: hbar           = 0.6582119514467406_dp     ! eV * fs
    real(dp), parameter :: twelfth        = 1.0_dp/12.0_dp
    integer,  parameter :: dimensionality = 3

    ! file units
    !integer, parameter  :: inp_unit = 67
    !integer, parameter  :: inp_unit = 38
    !integer, parameter  :: pes_unit = 38
    !integer, parameter  :: geo_unit = 55
    !integer, parameter  :: out_unit = 86
    !integer, parameter  :: inpnn_unit       = 61
    !integer, parameter  :: scaling_unit     = 62
    !integer, parameter  :: scalinge_unit     = 63
    !integer, parameter  :: weight_unit      = 64
    !integer, parameter  :: weighte_unit      = 65


    ! PES-related constants
    integer, parameter :: number_of_pess        = 8

    integer, parameter :: nparams_lj            = 2
    integer, parameter :: nparams_morse         = 3
    integer, parameter :: nparams_emt           = 7
    integer, parameter :: nparams_rebo          = 13
    integer, parameter :: nparams_nene          = 2

    integer, parameter :: pes_id_lj             = 2001
    integer, parameter :: pes_id_morse          = 2002
    integer, parameter :: pes_id_emt            = 2003
    integer, parameter :: pes_id_rebo           = 2004
    integer, parameter :: pes_id_ho             = 2005
    integer, parameter :: pes_id_simple_lj      = 2006
    integer, parameter :: pes_id_no_interaction = 2007
    integer, parameter :: pes_id_nene           = 2008
    integer, parameter :: pes_id_mace           = 2009

    ! PES names
    character(len=*), parameter :: pes_name_lj             = "lj"
    character(len=*), parameter :: pes_name_simple_lj      = "slj"
    character(len=*), parameter :: pes_name_emt            = "emt"
    character(len=*), parameter :: pes_name_rebo           = "rebo"
    character(len=*), parameter :: pes_name_no_interaction = "non"
    character(len=*), parameter :: pes_name_ho             = "ho"
    character(len=*), parameter :: pes_name_morse          = "morse"
    character(len=*), parameter :: pes_name_nene           = "nene"
    character(len=*), parameter :: pes_name_mace           = "mace"

    ! Debug ID
    integer, parameter :: debug_id_lj                 = 1
    integer, parameter :: debug_id_simple_lj          = 2
    integer, parameter :: debug_id_emt                = 3
    integer, parameter :: debug_id_rebo               = 4
    integer, parameter :: debug_id_no_interaction     = 5
    integer, parameter :: debug_id_ho                 = 6
    integer, parameter :: debug_id_morse              = 7
    integer, parameter :: debug_id_nene               = 8

     ! Internal program constants
    integer, parameter :: randseed(13)            = [7,5,3,11,9,1,17,2,9,6,4,5,8]
    integer, parameter :: max_string_length       = 1000
    real(dp), parameter :: tolerance              = 1.0e-9_dp


    ! Defaults
    integer,   parameter  :: default_int          = huge(0_4)
    real(dp),  parameter  :: default_real         = huge(0_dp)
    character, parameter  :: default_string       = ""
    logical,   parameter  :: default_bool         = .false.


    ! Propagation
    integer, parameter :: prop_id_verlet          = 1001
    integer, parameter :: prop_id_beeman          = 1002
    integer, parameter :: prop_id_langevin        = 1003
    integer, parameter :: prop_id_langevin_series = 1004
    integer, parameter :: prop_id_andersen        = 1005
    integer, parameter :: prop_id_pile            = 1006

    integer, parameter :: energy_and_force        = 3001
    integer, parameter :: energy_only             = 3002


    ! Geometry optimization
    integer, parameter :: geometry_opt_fire       = 4001


    ! Output
    character(len=*), parameter :: output_key_xyz         = "xyz"
    character(len=*), parameter :: output_key_energy      = "energy"
    character(len=*), parameter :: output_key_poscar      = "poscar"
    character(len=*), parameter :: output_key_vasp        = "vasp"
    character(len=*), parameter :: output_key_mxt         = "mxt"
    character(len=*), parameter :: output_key_scatter     = "scatter"
    character(len=*), parameter :: output_key_nene        = "nene"
    character(len=*), parameter :: output_key_aims        = "aims"
    character(len=*), parameter :: output_key_runner      = "runner"
    character(len=*), parameter :: output_key_is_adsorbed = "adsorption_status"
    character(len=*), parameter :: output_key_beads       = "beads"

    integer, parameter :: output_id_xyz         = 1
    integer, parameter :: output_id_energy      = 2
    integer, parameter :: output_id_poscar      = 3
    integer, parameter :: output_id_vasp        = 4
    integer, parameter :: output_id_mxt         = 5
    integer, parameter :: output_id_scatter     = 6
    integer, parameter :: output_id_nene        = 7
    integer, parameter :: output_id_aims        = 8
    integer, parameter :: output_id_runner      = 9
    integer, parameter :: output_id_is_adsorbed = 10
    integer, parameter :: output_id_beads       = 11

    ! Conversion constants to program units
    !
    ! Program basic units
    !           Length : Ang
    !           Time   : fs
    !           Energy : eV
    ! Program derived units
    !           Mass   : eV fs^2 / A^2 = 1/103.6382 amu
    !           Angle  : radian = 180 deg
    !           bohr   : bohr = 0.5291772 Angstroem
    real(dp), parameter :: amu2mass             = 103.638239276_dp
    real(dp), parameter :: deg2rad              = pi/180.0_dp
    real(dp), parameter :: rad2deg              = 180.0_dp/pi
    real(dp), parameter :: bohr2ang             = 0.529177211_dp
    real(dp), parameter :: ang2bohr             = 1.0_dp/bohr2ang
    real(dp), parameter :: p2GPa                = 160.2176565_dp
    real(dp), parameter :: joule2ev             = 6.2415093433e+18_dp
    real(dp), parameter :: kelvin2ev            = kB
    real(dp), parameter :: ha2ev                = 27.21138602_dp ! convert Ha to eV
    real(dp), parameter :: ev2ha                = 1.0_dp/ha2ev
    real(dp), parameter :: habohr2evang         = ha2ev*ang2bohr !51.4220670398_dp; convert Ha/bohr to eV/ang
    real(dp), parameter :: evang2habohr         = 1.0_dp/habohr2evang
    real(dp), parameter :: habohrcub2evangcub   = ha2ev*(ang2bohr**3)
    real(dp), parameter :: au2gpa               = 29419.844_dp ! 1 Ha/Bohr3 = 29419.844 GPa


    ! Atomic elements
    type(atomic_element), parameter :: elements(118) = [ &
                atomic_element('H  ',   1,      1.007), & ! Hydrogen
                atomic_element('He ',   2,      4.002), & ! Helium
                atomic_element('Li ',   3,      6.941), & ! Lithium
                atomic_element('Be ',   4,      9.012), & ! Beryllium
                atomic_element('B  ',   5,      10.811), & ! Boron
                atomic_element('C  ',   6,      12.011), & ! Carbon
                atomic_element('N  ',   7,      14.007), & ! Nitrogen
                atomic_element('O  ',   8,      15.999), & ! Oxygen
                atomic_element('F  ',   9,      18.998), & ! Fluorine
                atomic_element('Ne ',   10,     20.18), & ! Neon
                atomic_element('Na ',   11,     22.99), & ! Sodium
                atomic_element('Mg ',   12,     24.305), & ! Magnesium
                atomic_element('Al ',   13,     26.982), & ! Aluminum
                atomic_element('Si ',   14,     28.086), & ! Silicon
                atomic_element('P  ',   15,     30.974), & ! Phosphorus
                atomic_element('S  ',   16,     32.065), & ! Sulfur
                atomic_element('Cl ',   17,     35.453), & ! Chlorine
                atomic_element('Ar ',   18,     39.948), & ! Argon
                atomic_element('K  ',   19,     39.098), & ! Potassium
                atomic_element('Ca ',   20,     40.078), & ! Calcium
                atomic_element('Sc ',   21,     44.956), & ! Scandium
                atomic_element('Ti ',   22,     47.867), & ! Titanium
                atomic_element('V  ',   23,     50.942), & ! Vanadium
                atomic_element('Cr ',   24,     51.996), & ! Chromium
                atomic_element('Mn ',   25,     54.938), & ! Manganese
                atomic_element('Fe ',   26,     55.845), & ! Iron
                atomic_element('Co ',   27,     58.933), & ! Cobalt
                atomic_element('Ni ',   28,     58.693), & ! Nickel
                atomic_element('Cu ',   29,     63.546), & ! Copper
                atomic_element('Zn ',   30,     65.38), & ! Zinc
                atomic_element('Ga ',   31,     69.723), & ! Gallium
                atomic_element('Ge ',   32,     72.64), & ! Germanium
                atomic_element('As ',   33,     74.922), & ! Arsenic
                atomic_element('Se ',   34,     78.96), & ! Selenium
                atomic_element('Br ',   35,     79.904), & ! Bromine
                atomic_element('Kr ',   36,     83.798), & ! Krypton
                atomic_element('Rb ',   37,     85.468), & ! Rubidium
                atomic_element('Sr ',   38,     87.62), & ! Strontium
                atomic_element('Y  ',   39,     88.906), & ! Yttrium
                atomic_element('Zr ',   40,     91.224), & ! Zirconium
                atomic_element('Nb ',   41,     92.906), & ! Niobium
                atomic_element('Mo ',   42,     95.96), & ! Molybdenum
                atomic_element('Tc ',   43,     9855), & ! Technetium
                atomic_element('Ru ',   44,     101.07), & ! Ruthenium
                atomic_element('Rh ',   45,     102.906), & ! Rhodium
                atomic_element('Pd ',   46,     106.42), & ! Palladium
                atomic_element('Ag ',   47,     107.868), & ! Silver
                atomic_element('Cd ',   48,     112.411), & ! Cadmium
                atomic_element('In ',   49,     114.818), & ! Indium
                atomic_element('Sn ',   50,     118.71), & ! Tin
                atomic_element('Sb ',   51,     121.76), & ! Antimony
                atomic_element('Te ',   52,     127.6), & ! Tellurium
                atomic_element('I  ',   53,     126.904), & ! Iodine
                atomic_element('Xe ',   54,     131.293), & ! Xenon
                atomic_element('Cs ',   55,     132.905), & ! Cesium
                atomic_element('Ba ',   56,     137.327), & ! Barium
                atomic_element('La ',   57,     138.905), & ! Lanthanum
                atomic_element('Ce ',   58,     140.116), & ! Cerium
                atomic_element('Pr ',   59,     140.908), & ! Praseodymium
                atomic_element('Nd ',   60,     144.242), & ! Neodymium
                atomic_element('Pm ',   61,     145.00), & ! Promethium
                atomic_element('Sm ',   62,     150.36), & ! Samarium
                atomic_element('Eu ',   63,     151.964), & ! Europium
                atomic_element('Gd ',   64,     157.25), & ! Gadolinium
                atomic_element('Tb ',   65,     158.925), & ! Terbium
                atomic_element('Dy ',   66,     162.5), & ! Dysprosium
                atomic_element('Ho ',   67,     164.93), & ! Holmium
                atomic_element('Er ',   68,     167.259), & ! Erbium
                atomic_element('Tm ',   69,     168.934), & ! Thulium
                atomic_element('Yb ',   70,     173.054), & ! Ytterbium
                atomic_element('Lu ',   71,     174.967), & ! Lutetium
                atomic_element('Hf ',   72,     178.49), & ! Hafnium
                atomic_element('Ta ',   73,     180.948), & ! Tantalum
                atomic_element('W  ',   74,     183.84), & ! Wolfram
                atomic_element('Re ',   75,     186.207), & ! Rhenium
                atomic_element('Os ',   76,     190.23), & ! Osmium
                atomic_element('Ir ',   77,     192.217), & ! Iridium
                atomic_element('Pt ',   78,     195.084), & ! Platinum
                atomic_element('Au ',   79,     196.967), & ! Gold
                atomic_element('Hg ',   80,     200.59), & ! Mercury
                atomic_element('Tl ',   81,     204.383), & ! Thallium
                atomic_element('Pb ',   82,     207.20), & ! Lead
                atomic_element('Bi ',   83,     208.98), & ! Bismuth
                atomic_element('Po ',   84,     210.00), & ! Polonium
                atomic_element('At ',   85,     210.00), & ! Astatine
                atomic_element('Rn ',   86,     222.00), & ! Radon
                atomic_element('Fr ',   87,     223.00), & ! Francium
                atomic_element('Ra ',   88,     226.00), & ! Radium
                atomic_element('Ac ',   89,     227.00), & ! Actinium
                atomic_element('Th ',   90,     232.00), & ! Thorium
                atomic_element('Pa ',   91,     231.00), & ! Protactinium
                atomic_element('U  ',   92,     238.00), & ! Uranium
                atomic_element('Np ',   93,     237.00), & ! Neptunium
                atomic_element('Pu ',   94,     244.00), & ! Plutonium
                atomic_element('Am ',   95,     243.00), & ! Americium
                atomic_element('Cm ',   96,     247.00), & ! Curium
                atomic_element('Bk ',   97,     247.00), & ! Berkelium
                atomic_element('Cf ',   98,     251.00), & ! Californium
                atomic_element('Es ',   99,     252.00), & ! Einsteinium
                atomic_element('Fm ',   100,    257.00), & ! Fermium
                atomic_element('Md ',   101,    258.00), & ! Mendelevium
                atomic_element('No ',   102,    259.00), & ! Nobelium
                atomic_element('Lr ',   103,    262.00), & ! Lawrencium
                atomic_element('Rf ',   104,    261.00), & ! Rutherfordium
                atomic_element('Db ',   105,    262.00), & ! Dubnium
                atomic_element('Sg ',   106,    266.00), & ! Seaborgium
                atomic_element('Bh ',   107,    264.00), & ! Bohrium
                atomic_element('Hs ',   108,    267.00), & ! Hassium
                atomic_element('Mt ',   109,    268.00), & ! Meitnerium
                atomic_element('Ds ',   110,    271.00), & ! Darmstadtium
                atomic_element('Rg ',   111,    272.00), & ! Roentgenium
                atomic_element('Cn ',   112,    285.00), & ! Copernicium
                atomic_element('Nh ',   113,    284.00), & ! Nihonium
                atomic_element('Fl ',   114,    289.00), & ! Flerovium
                atomic_element('Mc ',   115,    288.00), & ! Moscovium
                atomic_element('Lv ',   116,    292.00), & ! Livermorium
                atomic_element('Ts ',   117,    295.00), & ! Tennessine
                atomic_element('Og ',   118,    294.00) & ! Oganesson
]

end module constants
