program hi
  use myKinds
  use pa
  use cmp
  use fft
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! declare variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: N = 4096
  integer     :: i,j
  integer     :: K = 0
  integer     :: Kneu = 100
  real(dp)    :: min = 0d0                                            ! minimum of energy in input file
  real(dp)    :: max = 0d0                                            ! maximum of energy in input file
  real(dp)    :: tau = 0d0                                            ! step width of fft-transformed signal
  real(dp)    :: step = 0d0                                           ! step width of input signal
  real(dp)    :: prec_deg = 1d-5                                      !
  integer     :: shift = 1                                            ! shift-parameter for cmp-module
  real(dp)    :: prec_cmp = 1d-5                                      ! precision-parameter for cmp_module
  complex(dp), dimension(2*N) :: c = (0d0,0d0)
  complex(dp), dimension(2*N) :: spektrum = (0d0,0d0)
  complex(dp), dimension(N)   :: a = (0d0,0d0)
  complex(dp), dimension(N)   :: b = (0d0,0d0)
  complex(dp), dimension(N)   :: z = (0d0,0d0)
  complex(dp), dimension(N)   :: d = (0d0,0d0)
  complex(dp), dimension(N)   :: w = (0d0,0d0)
  integer,     dimension(N)   :: zviel = 1

  integer :: io_error

  real(dp) :: readmax = 0d0                                           ! dummy array for input
  real(dp), dimension(2*N) :: readspektrum1 = 0d0                     ! dummy array for input
  real(dp), dimension(2*N) :: readspektrum2 = 0d0                     ! dummy array for input
  real(dp), dimension(2*N) :: array_out = 0d0                         ! dummy array for output

  character(len=64) :: in,sig,out                                     ! input and output files
  logical :: INPUTISREAL = .FALSE.                                    ! is input signal real?
  logical :: FOURIERAMPLIT = .FALSE.                                  ! should the amplitudes of the fft-transformed signal be given out?
  logical :: FOURIER = .TRUE.                                         ! should the signal be fft-transformed?
  logical :: PRINTSIGNAL = .FALSE.                                    ! should the fft-transformed signal be printed?
  logical :: MITTELWERT = .TRUE.                                      ! should the mean of the degenerated resonances be taken?
  logical :: GEOMETRIC = .FALSE.                                      ! should the mean be the geometric mean instead of the arithmetic?
  logical :: USE_CMP = .FALSE.                                        ! should the analysation be done twice and the results be compared?

  logical :: EXISTS                                                   ! for detecting the existence of files

  !Für die Kommandozeilenoptionen
  character(len=64) :: arg,arg2
  character(len=1)  :: testcharacter

  !Standardwerte setzen
  in = 'data.dat'
  sig = 'signal.dat'
  out = 'ergebnis.dat'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! parse command line options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  i = 0
  do while(i < command_argument_count())
!!$  do while(i < iargc())

     i = i + 1
     call get_command_argument(i, arg)
!!$     call getarg(i, arg)

     select case (arg)

     ! help
     case('-h', '--help')
        call print_help()
        call exit(0)

     ! real or complex input
     case('--real-input', '--no-complex-input')
        INPUTISREAL = .TRUE.
     case('--complex-input', '--no-real-input')
        INPUTISREAL = .FALSE.

     ! output file
     case('-o', '--output')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,out)
        read(out, "(1A1)") testcharacter
        if ( testcharacter == "-" ) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if

     ! input file
     case('-i', '--input')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,in)
        read(in, "(1A1)") testcharacter
        if ( testcharacter == "-" ) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if

     ! HI parameter
     case('-K')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,arg2)
        read(arg2, *, iostat=io_error) Kneu
        if(io_error /= 0) then
           write(0,*) "error: option ",trim(arg)," needs to be an integer"
           call exit(11)
        end if
        if(Kneu <= 0) then
           write(0,*) "error: argument of ",trim(arg)," needs to be an integer greater than zero"
           call exit(11)
        end if
        if(Kneu < 50) then
           write(0,*) "warning: argument of ",trim(arg)," is very small"
        end if
        if(Kneu > 200) then
           write(0,*) "warning: argument of ",trim(arg)," is very large"
        end if

     case('-g', '--degeneration-precision')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,arg2)
        read(arg2, *, iostat=io_error) prec_deg
        if(io_error /= 0) then
           write(0,*) "error: argument of option ",trim(arg)," needs to be a floating point number"
           call exit(11)
        end if
        if(prec_deg <= 0) then
           write(0,*) "error: argument of option ",trim(arg)," needs to be a floating point number greater than zero"
           call exit(11)
        end if
        if(prec_deg > 1d1) then
           write(0,*) "error: argument of option ",trim(arg)," is supposed to be too large"
           call exit(11)
        end if
        if(prec_deg > 1d-2) then
           write(0,*) "warning: argument of option ",trim(arg)," is very large"
        elseif(prec_deg < 1d-9) then
           write(0,*) "warning: argument of option ",trim(arg)," is very small"
        end if

     ! form of amplitudes
     case('--normal-amplitudes')
        FOURIERAMPLIT = .FALSE.
     case('--fourier-amplitudes')
        FOURIERAMPLIT = .TRUE.

     ! FFT?
     case('--no-fft', '--no-fourier', '--no-fourier-transformation')
        FOURIER = .FALSE.
        FOURIERAMPLIT = .TRUE.
     case('--fft', '--fourier', '--fourier-transformation')
        FOURIER = .TRUE.
        FOURIERAMPLIT = .FALSE.

     ! FFT-transformed signal?
     case('-s', '--signal')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,sig)
        read(sig, "(1A1)") testcharacter
        if ( testcharacter == "-" ) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        PRINTSIGNAL = .TRUE.
     case('--no-signal')
        PRINTSIGNAL = .FALSE.

     ! Mean
     case('--mean')
        MITTELWERT = .TRUE.
     case('--no-mean')
        MITTELWERT = .FALSE.
     case('--arithmetic-mean')
        MITTELWERT = .TRUE.
        GEOMETRIC = .FALSE.
     case('--geometric-mean')
        MITTELWERT = .TRUE.
        GEOMETRIC = .TRUE.

     ! Compare
     case('--no-c')
        USE_CMP = .FALSE.
     case('-c')
        USE_CMP = .TRUE.
     case('-c1', '--c-parameter-1')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,arg2)
        read(arg2, *, iostat=io_error) prec_cmp
        if(io_error /= 0) then
           write(0,*) "error: argument of option ",trim(arg)," needs to be a floating point number"
           call exit(11)
        end if
        if(prec_cmp < 0) then
           write(0,*) "error: argument of option ",trim(arg)," needs to be a non-negative floating point number"
           call exit(11)
        end if
!!$        if(prec_cmp > 1d-5) then
!!$           write(0,*) "warning: argument of option ",trim(arg)," is very large"
!!$        end if
        USE_CMP = .TRUE.
     case('-c2', '--c-parameter-2')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,arg2)
        read(arg2, *, iostat=io_error) shift
        if(io_error /= 0) then
           write(0,*) "error: option ",trim(arg)," needs to be an integer"
           call exit(11)
        end if
        if(shift < 0) then
           write(0,*) "error: argument of ",trim(arg)," needs to be an non-negative integer"
           call exit(11)
        end if
        USE_CMP = .TRUE.

     ! -- to end the command line argument parsing
     case('--')
        exit

     ! treat unknown command line argument
     case default
        write(0,*) "error: unknown command line option ",trim(arg)
        call exit(11)

     end select

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! treat command line arguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! --normal-amplitudes and --no-fft should not be given together
  if( (.not. FOURIERAMPLIT) .and. (.not. FOURIER)) then
     write(0,*) "error: conflicting command line arguments:"
     write(0,*) "       normal amplitudes can't be computed in no-fft mode"
     call exit(1)
  end if
  ! -s and --no-fft should not be given together
  if(PRINTSIGNAL .and. (.not. FOURIER)) then
     write(0,*) "error: conflicting command line arguments:"
     write(0,*) "       no fft-transformation should be done but fft-transformed signal should be printed"
     call exit(1)
  end if
  ! parameter for cmp such that cmp doesn't do anything
  if(prec_cmp == 0d0) then
     write(0,*) "warning: argument of option -c1 resp. --c-parameter-1 is zero; ignoring -c"
     USE_CMP = .FALSE.
  end if
  if(shift == 0) then
     write(0,*) "warning: argument of option -c2 resp. --c-parameter-2 is zero; ignoring -c"
     USE_CMP = .FALSE.
  end if
  ! do the files exist?
  inquire(file=in, exist=EXISTS)
  if(.not. EXISTS) then
     write(0,*) "error: input file '",trim(in),"' does not exist"
     call exit(1)
  end if
  inquire(file=out, exist=EXISTS)
  if(EXISTS) then
     write(0,*) "error: output file '",trim(out),"' does already exist"
     call exit(1)
  end if
  if(PRINTSIGNAL) then
     inquire(file=sig, exist=EXISTS)
     if(EXISTS) then
        write(0,*) "error: output file '",trim(sig),"' does already exist"
        call exit(1)
     end if
  end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! read the signal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(unit=11, file=in, status='old', action='read', iostat=io_error)
  if(io_error /= 0) then
     write(*,*) "error",io_error,"while opening the file '",trim(in),"'"
     call exit(1)
  end if
  if(INPUTISREAL) then
     ! assuming the input signal is real
     K = 1
     read(11, *, iostat=io_error) min, readspektrum1(K)
     ! treat error
     if(io_error /= 0) then
        write(0,*) "error: no data in '",trim(in),"'"
        call exit(1)
     end if
     ! read rest of the data
     max = min
     do
        if(K+1 > 2*N) then
           write(*,*) "warning: too many data points, using just the first",K
           exit
        end if
        read(11, *, iostat=io_error) readmax, readspektrum1(K+1)
        if(io_error > 0) then
           write(*,*) "error",io_error,"in line",K+1,"while reading the input file"
           exit
        end if
        K = K + 1
        max = readmax
     end do
  else
     ! assuming the input signal is complex
     K = 1
     read(11, *, iostat=io_error) min, readspektrum1(K), readspektrum2(K)
     ! treat error
     if(io_error /= 0) then
        write(0,*) "error: no data in '",trim(in),"'"
        call exit(1)
     end if
     ! read rest of the data
     max = min
     do
        if(K+1 > 2*N) then
           write(*,*) "warning: too many data points, using just the first",K
           exit
        end if
        read(11, *, iostat=io_error) readmax, readspektrum1(K+1), readspektrum2(K+1)
        if(io_error < 0) then
           exit
        elseif(io_error > 0) then
           write(*,*) "error",io_error,"in line",K+1,"while reading the input file"
           exit
        end if
        K = K + 1
        max = readmax
     end do
  end if
  close(11)

  forall(i = 1:K) spektrum(i) = cmplx(readspektrum1(i),readspektrum2(i))
  step = (max-min)/(K-1)
  tau = step
  max = max + step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! treat illegal values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! max > min
  if(step <= 0d0) then
     write(0,*) "error: max <= min"
     call exit(11)
  end if
  ! K not to small
  if(K < 128) then
     write(*,*) "warning: the input signal is very short (just",K,"data points)"
  end if
  ! in no-fft mode: min = 0
  ! Verbessern? abs(min) > ???
  if( (abs(min) > epsilon(max)) .and. (.not. FOURIER) ) then
     write(*,*) "error: min ist not zero in no-fft mode; program must be extended"
     call exit(11)
  end if
  ! K a power of 2 because of FFT
  if(FOURIER) then
     ! K > 0 after reading the signal
     if( IAND(K,(K-1)) /= 0 ) then
        K = 2**(int( floor( log(dble(K))/log(2d0) ) ))
        max = min + step*K
        write(*,*) "warning: the number of signal points of the input signal is no power of 2"
        write(*,*) "         using only the first",K
     end if
  end if
  ! Kneu <= K/2
  if(Kneu > K/2) then
     Kneu = K/2
     write(*,*) "warning: argument of -K is too large; using",Kneu,"instead"
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! harmonic-inversion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! fft-part
  if(FOURIER) then
     ! fft-transformation
     call fft_fourier(min,max,N,K,spektrum,c,tau)
     ! print fft-transformed signal
     if(PRINTSIGNAL) then
        open(unit=31, file=sig, status='replace', action='write', iostat=io_error)
        if(io_error /= 0) then
           write(0,*) "error", io_error, "while creating the file '",trim(sig),"'"
           call exit(1)
        end if
        do i=1, K
           write(31,*) c(i)
        end do
        close(31)
     end if
  else
     c(1:K) = spektrum(1:K)
  end if

  !Signal analysieren und Energien und Amplituden auslesen
  if(USE_CMP) then
     call cmp_all(N,K,Kneu,shift,prec_cmp,tau,prec_deg,MITTELWERT,GEOMETRIC,c,z,zviel,w,d)
  else
     call pa_all(N,Kneu,tau,prec_deg,MITTELWERT,GEOMETRIC,c,z,zviel,w,d)
  end if
  
  !Amplituden und Energien richtig skalieren
  if(FOURIERAMPLIT) then
     call pa_skal_fourier(N,Kneu,min,max,tau,zviel,w,d)
  else
     call pa_skal(N,Kneu,min,max,tau,zviel,w,d)
  end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Verbessern:
!!! Ausgabe?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Ergebnisse ausgeben
  open(unit=21, file=out, status='replace', action='write', iostat=io_error)
  if(io_error /= 0) then
     write(0,*) "error",io_error,"while creating the file '",trim(out),"'"
     call exit(1)
  end if

  i = 1
  do while (i <= Kneu)
     if(zviel(i) <= 0) then
        i = i + 1
     else
        do j=0,zviel(i)-1
           array_out(2*zviel(i)-1-2*j) = dble(d(i+j))
           array_out(2*zviel(i)-2*j) = dimag(d(i+j))
        end do
        write(21,*) zviel(i),dble(w(i)),dimag(w(i)),array_out(1:2*zviel(i))
        i = i + zviel(i)
     end if
  end do

  close(21)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  ! print help message
  subroutine print_help()
    implicit none

    write(*,*) "PROGRAM: harmonic-inversion"
    write(*,*) "VERSION: 14 (25.03.2014)"
    write(*,*) "this program analyses a resonance spektra with the harmonic-inversion-method"
    write(*,*) ""
    write(*,*) "OPTIONS:"
    write(*,*) "-h, --help"
    write(*,*) "        print this help and exit"
    write(*,*) "-i <file>, --input <file>"
    write(*,*) "        set the input file to <file>"
    write(*,*) "        the file must contain the energy (first column) and the spectra (next column(s))"
    write(*,*) "        if not referred, 'data.dat' ist used"
    write(*,*) "--real-input, --no-complex-input"
    write(*,*) "        assumes real input signal"
    write(*,*) "        overrides previous --complex-input options"
    write(*,*) "--complex-input, --no-real-input"
    write(*,*) "        assumes complex input signal"
    write(*,*) "        overrides previous --real-input options"
    write(*,*) "        enabled by default"
    write(*,*) "-o <file>, --output <file>"
    write(*,*) "        set the output file to <file>"
    write(*,*) "        it contains results of the calculation"
    write(*,*) "        if not referred, 'ergebnis.dat' ist used"
    write(*,*) "--normal-amplitudes, --no-fourier-amplitudes"
    write(*,*) "        normal form of amplitudes are given out"
    write(*,*) "        overrides a previous --fourier-amplitudes options"
    write(*,*) "--fourier-amplitudes"
    write(*,*) "        amplitudes are given out as amplitudes of the fourier-transformed signal"
    write(*,*) "        overrides previous --normal-amplitudes options"
    write(*,*) "--no-fft, --no-fourier, --no-fourier-transformation"
    write(*,*) "        disable fourier-transformation of the input signal"
    write(*,*) "        overrides previous --fft options"
    write(*,*) "        implys the --fourier-amplitudes option"
    write(*,*) "--fft, --fourier, --fourier-transformation"
    write(*,*) "        enable fourier-transformation of the input signal"
    write(*,*) "        overrides a previous --no-fft options"
    write(*,*) "        implys the --normal-amplitudes options"
    write(*,*) "        enabled by default"
    write(*,*) "--no-signal"
    write(*,*) "        the fourier-transformed data is not given out"
    write(*,*) "        overrides a previous --signal option"
    write(*,*) "--signal"
    write(*,*) "        the fourier-transformed data is given out"
    write(*,*) "        overrides a previous --no-signal options"
    write(*,*) "        enabled by default"
    write(*,*) "-s <file>, --signal-file <file>"
    write(*,*) "        set output file of the fourier-transformed data to <file>"
    write(*,*) "        if not referred, 'signal.dat' ist used"
    write(*,*) "        implys the --signal option"
    write(*,*) "-K <integer>"
    write(*,*) "        the number of points which are used to calculate the energies and amplitudes"
    write(*,*) "        <integer> should be an integer between 50 and 200,"
    write(*,*) "            but it mustn't be larger than the number of input data points"
    write(*,*) "        if not referred, 100 is used"
    write(*,*) "-g <float>, --degeneration-precision <float>"
    write(*,*) "        the precision which is used to identify exceptional points"
    write(*,*) "        <float> should be a floating point number between 1e-3 and 1e-7,"
    write(*,*) "            but it mustn't be greater than zero"
    write(*,*) "        if not referred, 1e-5 is used"
    write(*,*) "--mean"
    write(*,*) "        take the mean of the degenerated resonances"
    write(*,*) "        overrides a previous --no-mean option"
    write(*,*) "        enabled by default"
    write(*,*) "--no-mean"
    write(*,*) "        take no mean of the degenerated resonances"
    write(*,*) "        overrides a previous --mean option"
    write(*,*) "--arithmetic-mean"
    write(*,*) "        compute the arithmetic mean of the z-parameter of the degenerated resonances"
    write(*,*) "        this option implys the --mean option"
    write(*,*) "        implys the --mean option"
    write(*,*) "        enabled by default"
    write(*,*) "--geometric-mean"
    write(*,*) "        compute the geometric mean of the z-parameter of the degenerated resonances"
    write(*,*) "            (this is equivalent of computing the arithmetic mean of the w-parameters)"
    write(*,*) "        implys the --mean option"
    write(*,*) "--no-c"
    write(*,*) "        disables the -c option"
    write(*,*) "        overrides previous -c option"
    write(*,*) "        enabled by default"
    write(*,*) "-c"
    write(*,*) "        the results are computed twice - one time with beginning with the first signal point, one time not"
    write(*,*) "        the two results are compared and just the points occuring in both results are given out"
    write(*,*) "        overrides a previous --no-c option"
    write(*,*) "-c1 <float>, --c-parameter-1 <float>"
    write(*,*) "        specifies the maximum relative difference of the points in both results in terms of machine epsilon"
    write(*,*) "            (when -c option is enabled)"
    write(*,*) "        implys the -c option"
    write(*,*) "        <float> should be a non-negative floating point number"
    write(*,*) "        if not referred, ten is used"
    write(*,*) "-c2 <integer>, --c-parameter-2 <integer>"
    write(*,*) "        specifies the first signal point (start counting by 0) used for the second computation"
    write(*,*) "            (when -c option is enabled)"
    write(*,*) "        implys the -c option"
    write(*,*) "        <integer> should be a non-negative integer"
    write(*,*) "        if not referred, 1 is used"

  end subroutine print_help

end program hi
