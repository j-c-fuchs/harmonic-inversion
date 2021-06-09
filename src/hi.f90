! harmonic inversion
! Program for harmonic inversion analysis
!
! Copyright (C) 2014  Jacob Fuchs
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

program hi
  use myKinds
  use pa
  use fft
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! main program for the harmonic inversion analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! some general comments:
!!!   - the program offers command line arguments; the command line
!!!     argument -h will show a short help message
!!!   - the important routines for the harmonic inversion analysis are
!!!     in the modules fft and pa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! use static arrays -- N is the length of the static arrays
  integer, parameter :: N = 4096
  integer     :: i,j
  integer     :: len = 0           ! signal length
  integer     :: K = 100           ! number of resonances to extract
  real(dp)    :: min = 0_dp        ! minimum value of input signal
  real(dp)    :: max = 0_dp        ! maximum value of input signal
  real(dp)    :: tau = 0_dp        ! time step tau
  real(dp)    :: step = 0_dp       ! step width of input signal
  ! maximum distance between to resonances to treat them as degenerated
  real(dp)    :: prec_deg = 1d-5
  ! input signal or spectra in the energy-domain, respectively
  complex(dp), dimension(2*N) :: spectra = (0_dp,0_dp)
  ! signal to analyse
  complex(dp), dimension(2*N) :: c = (0_dp,0_dp)
  ! wanted z's, w's and d's
  complex(dp), dimension(N)   :: z = (0_dp,0_dp)
  complex(dp), dimension(N)   :: d = (0_dp,0_dp)
  complex(dp), dimension(N)   :: w = (0_dp,0_dp)
  ! corresponding degree of degeneration
  integer,     dimension(N)   :: z_dod = 1

  integer :: io_error

  real(dp) :: readmax = 0_dp       ! dummy variable for input
  ! dummy arrays for input
  real(dp), dimension(2*N) :: readspectra1 = 0_dp
  real(dp), dimension(2*N) :: readspectra2 = 0_dp
  ! dummy array for output
  real(dp), dimension(2*N) :: array_out = 0_dp

  ! input and output files
  character(len=64) :: in,sig,out
  ! is input signal real (instead of complex)?
  logical :: INPUTISREAL = .FALSE.
  ! should the amplitudes be given out as coefficients of the
  ! fft-transformed signal (i.e. d~ in the paper, see Eqn. (36))
  logical :: FOURIERAMPLIT = .FALSE.
  ! should the signal be fft-transformed?
  logical :: FOURIER = .TRUE.
  ! should the fft-transformed signal be printed?
  logical :: PRINTSIGNAL = .FALSE.
  ! compute mean of the frequencies?
  ! (compute means of the z's of degenerated resonances)
  logical :: MEAN = .TRUE.
  ! compute geometric mean instead of the arithmetic of the
  ! frequencies?
  ! (taking the geometric mean means taking the arithmetic mean of the
  ! corresponding w's)
  logical :: GEOMETRIC = .FALSE.

  ! dummy variable for detecting the existence of files
  logical :: EXISTS

  ! for parsing the command line arguments
  character(len=64) :: arg,arg2
  character(len=1)  :: testcharacter

  ! set standard values for file names
  in = 'data.dat'
  sig = 'signal.dat'
  out = 'ergebnis.dat'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! parse command line options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        call get_command_argument(i, out)
!!$        call getarg(i,out)
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
        call get_command_argument(i, in)
!!$        call getarg(i,in)
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
        call get_command_argument(i, arg2)
!!$        call getarg(i,arg2)
        read(arg2, *, iostat=io_error) K
        if(io_error /= 0) then
           write(0,*) "error: option ",trim(arg)," needs to be an&
                & integer"
           call exit(11)
        end if
        if(K <= 0) then
           write(0,*) "error: argument of ",trim(arg)," needs to be an&
                & integer greater than zero"
           call exit(11)
        end if

     case('-g', '--degeneration-precision')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call get_command_argument(i, arg2)
!!$        call getarg(i,arg2)
        read(arg2, *, iostat=io_error) prec_deg
        if(io_error /= 0) then
           write(0,*) "error: argument of option ",trim(arg)," needs to&
                & be a floating point number"
           call exit(11)
        end if
        if(prec_deg <= 0) then
           write(0,*) "error: argument of option ",trim(arg)," needs to&
                & be a floating point number greater than zero"
           call exit(11)
        end if

     ! form of amplitudes to give out
     case('--normal-amplitudes')
        FOURIERAMPLIT = .FALSE.
     case('--fourier-amplitudes')
        FOURIERAMPLIT = .TRUE.

     ! perform the FFT?
     case('--no-fft', '--no-fourier', '--no-fourier-transformation')
        FOURIER = .FALSE.
        FOURIERAMPLIT = .TRUE.
     case('--fft', '--fourier', '--fourier-transformation')
        FOURIER = .TRUE.
        FOURIERAMPLIT = .FALSE.

     ! print out the FFT-transformed signal?
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
        MEAN = .TRUE.
     case('--no-mean')
        MEAN = .FALSE.
     case('--arithmetic-mean')
        MEAN = .TRUE.
        GEOMETRIC = .FALSE.
     case('--geometric-mean')
        MEAN = .TRUE.
        GEOMETRIC = .TRUE.

     ! -- to end the command line argument parsing
     case('--')
        exit

     ! treat unknown command line argument
     case default
        write(0,*) "error: unknown command line option ",trim(arg)
        call exit(11)

     end select

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! check for conflicting command line arguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! --normal-amplitudes and --no-fft should not be given together
  if( (.not. FOURIERAMPLIT) .and. (.not. FOURIER)) then
     write(0,*) "error: conflicting command line arguments:"
     write(0,*) "       normal amplitudes can't be computed in no-fft&
          & mode"
     call exit(1)
  end if

  ! -s and --no-fft should not be given together
  if(PRINTSIGNAL .and. (.not. FOURIER)) then
     write(0,*) "error: conflicting command line arguments:"
     write(0,*) "       no fft-transformation should be done but&
          & fft-transformed signal should be printed"
     call exit(1)
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
        write(0,*) "error: output file '",trim(sig),"' does already&
             & exist"
        call exit(1)
     end if
  end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! read the signal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! open input file
  open(unit=11, file=in, status='old', action='read', iostat=io_error)
  ! treat errors while opening the input file
  if(io_error /= 0) then
     write(*,*) "error",io_error,"while opening the file '",trim(in),"'"
     call exit(1)
  end if

  if(INPUTISREAL) then
     ! assuming the input signal is real
     len = 1
     ! read the first line
     read(11, *, iostat=io_error) min, readspectra1(len)
     ! treat errors while reading input
     if(io_error /= 0) then
        write(0,*) "error: no data in '",trim(in),"'"
        call exit(1)
     end if
     ! read rest of the data
     max = min
     do
        ! warning when more datapoints are available than necessary
        if(len+1 > 2*N) then
           write(*,*) "warning: too many data points, using just the&
                & first",len
           exit
        end if
        read(11, *, iostat=io_error) readmax, readspectra1(len+1)
        ! check if EOF is reached
        if(io_error < 0) exit
        ! treat errors
        if(io_error > 0) then
           write(*,*) "error",io_error,"in line",len+1,"while reading&
                & the input file"
           exit
        end if
        len = len + 1
        max = readmax
     end do
  else
     ! assuming the input signal is complex
     len = 1
     ! read the first line
     read(11, *, iostat=io_error) min, readspectra1(len),&
          & readspectra2(len)
     ! treat error
     if(io_error /= 0) then
        write(0,*) "error: no data in '",trim(in),"'"
        call exit(1)
     end if
     ! read rest of the data
     max = min
     do
        ! warning when more datapoints are available than necessary
        if(len+1 > 2*N) then
           write(*,*) "warning: too many data points, using just the&
                & first",len
           exit
        end if
        read(11, *, iostat=io_error) readmax, readspectra1(len+1),&
             & readspectra2(len+1)
        ! check if EOF is reached
        if(io_error < 0) exit
        ! treat errors
        if(io_error > 0) then
           write(*,*) "error",io_error,"in line",len+1,"while reading&
                & the input file"
           exit
        end if
        len = len + 1
        max = readmax
     end do
  end if

  ! close input file
  close(11)

  ! assign read signal to array spectra and compute step width etc.
  forall(i = 1:len)&
       & spectra(i) = cmplx(readspectra1(i), readspectra2(i), dp)
  step = (max-min)/(len-1)
  tau = step
  max = max + step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! treat illegal values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! maximum value <= minimum value
  if(step <= 0_dp) then
     write(0,*) "error: max <= min"
     call exit(11)
  end if

  ! in no-fft mode: min == 0
  if( (abs(min) > epsilon(max)) .and. (.not. FOURIER) ) then
     write(*,*) "error: min ist not zero in no-fft mode; program must&
          & be extended"
     call exit(11)
  end if

  ! is the signal length a power of 2? (this is better for the FFT)
  if(FOURIER) then
     ! len > 0 after reading the signal
     if( IAND(len,(len-1)) /= 0 ) then
        len = 2**(int( floor( log(dble(len))/log(2d0) ) ))
        max = min + step*len
        write(*,*) "warning: the number of signal points of the input&
             & signal is no power of 2"
        write(*,*) "         using only the first",len
     end if
  end if

  ! K requires at least twice as much signal points than existing
  if(K > len/2) then
     K = len/2
     write(*,*) "warning: argument of -K is too large;&
          & using",K,"instead"
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(FOURIER) then
     ! fft-transformation
     call fft_fourier(min,max,N,len,spectra,c,tau)
     ! print fft-transformed signal
     if(PRINTSIGNAL) then
        open(unit=31, file=sig, status='replace', action='write',&
             & iostat=io_error)
        if(io_error /= 0) then
           write(0,*) "error", io_error, "while creating the file&
                & '",trim(sig),"'"
           call exit(1)
        end if
        do i=1, len
           write(31,*) c(i)
        end do
        close(31)
     end if
  else
     c(1:len) = spectra(1:len)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! solve the equations using the Pade-approximation and bring the
!!! results into the right form
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! analysis
  call pa_all(N, K, tau, prec_deg, MEAN, GEOMETRIC,&
       & c, z, z_dod, w, d)
  
  ! scale the amplitudes and frequencies
  if(FOURIERAMPLIT) then
     call pa_skal_fourier(N, K, min, max, tau, z_dod, w, d)
  else
     call pa_skal(N, K, min, max, tau, z_dod, w, d)
  end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! open the output file
  open(unit=21, file=out, status='replace', action='write',&
       & iostat=io_error)
  ! checking for errors
  if(io_error /= 0) then
     write(0,*) "error",io_error,"while creating the file&
          & '",trim(out),"'"
     call exit(1)
  end if

  ! write results to output file
  i = 1
  do while (i <= K)
     ! avoid endless loops
     if(z_dod(i) <= 0) then
        i = i + 1
     else
        do j=0,z_dod(i)-1
           array_out(2*z_dod(i)-1-2*j) = real(d(i+j), dp)
           array_out(2*z_dod(i)-2*j) = dimag(d(i+j))
        end do
        write(21,*) z_dod(i), real(w(i), dp), dimag(w(i)),&
             & array_out(1:2*z_dod(i))
        i = i + z_dod(i)
     end if
  end do

  ! close the output file
  close(21)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine to print the help message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  ! print help message
  subroutine print_help()
    implicit none

    write(*,*) "PROGRAM: harmonic-inversion"
    write(*,*) "VERSION: revised 28.08.2014"
    write(*,*) "DESCRIPTION:"
    write(*,*) "this program analyses a resonance spektra with the&
         & harmonic-inversion-method"
    write(*,*) ""
    write(*,*) "INPUT FORMAT:"
    write(*,*) "the input file should contain following columns in&
         & following order:"
    write(*,*) "  (1) value of energy/time"
    write(*,*) "  (2) real part of signal"
    write(*,*) "  (3) (if necessary) imaginary part of signal"
    write(*,*) ""
    write(*,*) "OUTPUT FORMAT:"
    write(*,*) "the output file should contain following columns in&
         & following order:"
    write(*,*) "  (1) degree of degeneration"
    write(*,*) "  (2) real part of the resonance"
    write(*,*) "  (3) imaginary part of the resonance"
    write(*,*) "  (4) real part of first amplitude"
    write(*,*) "  (5) imaginary part of first amplitude"
    write(*,*) "  (6) real part of second amplitude (if degree of&
         & degeneration is higher than 1)"
    write(*,*) "  (7) imaginary part of second amplitude (if degree of&
         & degeneration is higher than 1)"
    write(*,*) "        further amplitudes if degree of degeneration is&
         & higher"
    write(*,*) ""
    write(*,*) "AVAILABLE COMMAND LINE ARGUMENTS:"
    write(*,*) "-h, --help"
    write(*,*) "        print this help and exit"
    write(*,*) "-i <file>, --input <file>"
    write(*,*) "        set the input file to <file>"
    write(*,*) "        the file must contain the energy (first column)&
         & and the spectra (next column(s))"
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
    write(*,*) "        overrides a previous --fourier-amplitudes&
         & options"
    write(*,*) "--fourier-amplitudes"
    write(*,*) "        amplitudes are given out as amplitudes of the&
         & fourier-transformed signal"
    write(*,*) "        overrides previous --normal-amplitudes options"
    write(*,*) "--no-fft, --no-fourier, --no-fourier-transformation"
    write(*,*) "        disable fourier-transformation of the input&
         & signal"
    write(*,*) "        overrides previous --fft options"
    write(*,*) "        implys the --fourier-amplitudes option"
    write(*,*) "--fft, --fourier, --fourier-transformation"
    write(*,*) "        enable fourier-transformation of the input&
         & signal"
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
    write(*,*) "        set output file of the fourier-transformed data&
         & to <file>"
    write(*,*) "        if not referred, 'signal.dat' ist used"
    write(*,*) "        implys the --signal option"
    write(*,*) "-K <integer>"
    write(*,*) "        the number of points which are used to&
         & calculate the energies and amplitudes"
    write(*,*) "        <integer> should be an integer between 50 and&
         & 200,"
    write(*,*) "            but it mustn't be larger than the number of&
         & input data points"
    write(*,*) "        if not referred, 100 is used"
    write(*,*) "-g <float>, --degeneration-precision <float>"
    write(*,*) "        the precision which is used to identify&
         & exceptional points"
    write(*,*) "        <float> should be a floating point number&
         & between 1e-3 and 1e-7,"
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
    write(*,*) "        compute the arithmetic mean of the z-parameter&
         & of the degenerated resonances"
    write(*,*) "        this option implys the --mean option"
    write(*,*) "        implys the --mean option"
    write(*,*) "        enabled by default"
    write(*,*) "--geometric-mean"
    write(*,*) "        compute the geometric mean of the z-parameter&
         & of the degenerated resonances"
    write(*,*) "            (this is equivalent of computing the&
         & arithmetic mean of the w-parameters)"
    write(*,*) "        implys the --mean option"

  end subroutine print_help

end program hi
