program hi
  use myKinds
  use pa
  use fft
  implicit none

  integer, parameter :: N = 4096
  integer     :: K,i,j,Kneu
  real(dp)    :: tau                                                  !Schrittweite im Zeitsignal
  real(dp)    :: min,max                                              !Minimale und Maximale Energie des Datensatzes
  real(dp)    :: schritt                                              !Schrittweite im Energiesignal
  real(dp)    :: genau                                                !Schrittweite im Energiesignal
  complex(dp), dimension(2*N) :: c,spektrum
  complex(dp), dimension(N)   :: a,b,z,d,w
  integer,     dimension(N)   :: zviel = 1

  integer :: io_error

  character(len=64) :: in1,in2,sig,out                                !In- und Output-Dateien
  logical :: PRINTSIGNAL = .TRUE.                                     !Soll das Fouriertransformierte Signal ausgegeben werden?

  !Für die Kommandozeilenoptionen
  character(len=64) :: arg,arg2
  character(len=1)  :: testcharacter

  !Standardwerte setzen
  Kneu = 100
  genau = 1d-5

  in1 = 'data.dat'
  in2 = 'data.dat2'
  sig = 'signal.dat'
  out = 'ergebnis.dat'

  !Parse command line options
  i = 0
  do while(i < iargc())
     i = i + 1
     call getarg(i, arg)
     select case (arg)
     case('-h', '--help')
        call print_help()
        call exit(0)
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
     case('-i1', '--input1')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,in1)
        read(in1, "(1A1)") testcharacter
        if ( testcharacter == "-" ) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
     case('-i2', '--input2')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,in2)
        read(in2, "(1A1)") testcharacter
        if ( testcharacter == "-" ) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
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
     case('--no-signal')
        PRINTSIGNAL = .FALSE.
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
     case('-g', '--genau')
        i = i + 1
        if ( i > iargc()) then
           write(0,*) "error: option ",trim(arg)," requires an argment"
           call exit(11)
        end if
        call getarg(i,arg2)
        read(arg2, *, iostat=io_error) genau
        if(io_error /= 0) then
           write(0,*) "error: argument of option ",trim(arg)," needs to be a floating point number"
           call exit(11)
        end if
        if(genau <= 0) then
           write(0,*) "error: argument of option ",trim(arg)," needs to be a floating point number greater than zero"
           call exit(11)
        end if
        if(genau > 1d-2) then
           write(0,*) "warning: argument of option ",trim(arg)," is very large"
        end if
        if(genau > 1d1) then
           write(0,*) "error: argument of option ",trim(arg)," is supposed to be too large"
           call exit(11)
        end if
        if(genau < 1d-9) then
           write(0,*) "warning: argument of option ",trim(arg)," is very small"
        end if
     case('--')
        exit
     case default
        write(0,*) "error: unknown command line option"
        call exit(11)
     end select
  end do

  !Read K and tau
  open(unit=12, file=in2, status='old', action='read', iostat=io_error)
  if(io_error /= 0) then
     write(0,*) "error", io_error, "while opening the file '",trim(in2),"'"
     call exit(1)
  end if
  read(12, *, iostat=io_error) K
  read(12, *, iostat=io_error) min
  read(12, *, iostat=io_error) max
  close(12)

  !Behandlung der Fehler und der unzulässigen Werte
  if(io_error /= 0) then
     write(0,*) "error while reading the file '",trim(in2),"'"
     call exit(1)
  end if
  if(K > N) then
     write(0,*) "error: too many data points"
     call exit(2)
  end if
  if(max <= min) then
     write(0,*) "error in data in file ",in2,": max <= min"
     call exit(1)
  end if
  if(Kneu > K) then
     write(0,*) "warning: argument of -K is too large"
     call exit(11)
  end if

  !Berechnen der Schrittweite
  schritt = (max-min)/K

  !Read the Signal
  open(unit=11, file=in1, status='old', action='read', iostat=io_error)
  if(io_error /= 0) then
     write(0,*) "error", io_error, "while opening the file '",trim(in1),"'"
     call exit(1)
  end if
  do i=1, K
     read(11, *) spektrum(i)
  end do
  close(11)

  !Signal Fouriertransformieren
  call fft_fourier(min,max,N,K,spektrum,c,tau)

  !Signal ausgeben, wenn gewünscht
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

  !Signal analysieren und Energien und Amplituden auslesen
  call pa_all(N,Kneu,tau,genau,c,z,zviel,w,d)

  !Amplituden und Energien richtig skalieren
  call pa_skal(N,Kneu,min,max,tau,zviel,w,d)

!!$  !Amplituden richtig skalieren
!!$  i = 1
!!$  do while(i <= Kneu)
!!$     !einfache Nullstelle
!!$     if(zviel(i) == 1) then
!!$        d(i) = (0d0,1d0) * d(i)
!!$     end if
!!$     !zweifache Nullstelle
!!$     if(zviel(i) == 2) then
!!$        d(i) = (-1d0/tau) * d(i)
!!$        d(i+1) = (0d0,1d0) * d(i+1)
!!$     end if
!!$     !höhere Nullstellen
!!$     if(zviel(i) > 2) then
!!$        do j=i, i+zviel(i)-1
!!$           d(j) = (0d0,0d0)
!!$        end do
!!$        write(*,*) "warning: order of degeneration higher than 2 for z = ", z(i)
!!$     end if
!!$     !Sicherheitsabfrage zur Vermeidung von Endlosschleifen
!!$     if(zviel(i) < 1) then
!!$        write(0,*) "error: Vielfachheit kleiner 1 (Berechnung der Amplituden)"
!!$        call exit(2)
!!$     end if
!!$     i = i + zviel(i)
!!$  end do
!!$
!!$  !Energien verschieben
!!$  do i=1, Kneu
!!$     if(dble(w(i)) < 0) w(i) = w(i) + 8d0*atan(1d0)
!!$     w(i) = w(i) + min
!!$  end do

  !Ergebnisse ausgeben
  open(unit=21, file=out, status='replace', action='write', iostat=io_error)
  if(io_error /= 0) then
     write(0,*) "error", io_error, "while creating the file '",trim(out),"'"
     call exit(1)
  end if

  i = 1
  do while (i <= Kneu)
     write(21,'(X,I3$)') zviel(i)
     write(21,*) w(i)
     do j=i, i+zviel(i)-1
        write(21,*) "   ", d(j)
     end do
     if(zviel(i) <= 0) then
        i = i + 1
     else
        i = i + zviel(i)
     end if
  end do

  close(21)

contains

  subroutine print_help()
    implicit none

    write(*,*) "PROGRAM: harmonic-inversion"
    write(*,*) "VERSION: 13 (12.07.2013)"
    write(*,*) "this program analyses a resonance spektra with the harmonic-inversion-method"
    write(*,*) ""
    write(*,*) "OPTIONS:"
    write(*,*) "-h, --help"
    write(*,*) "        print this help and exit"
    write(*,*) "-i1 <file>, --input1 <file>"
    write(*,*) "        set the first input file to <file>"
    write(*,*) "        the file must contain the spectra"
    write(*,*) "        if not referred, 'data.dat' ist used"
    write(*,*) "-i2 <file>, --input2 <file>"
    write(*,*) "        set the second input file to <file>"
    write(*,*) "        the file must contain the number of data points and the minimum and the maximum of the energy"
    write(*,*) "        if not referred, 'data.dat2' ist used"
    write(*,*) "-o <file>, --output <file>"
    write(*,*) "        set the output file to <file>"
    write(*,*) "        it contains results of the calculation"
    write(*,*) "        if not referred, 'ergebnis.dat' ist used"
    write(*,*) "-s <file>, --signal <file>"
    write(*,*) "        set output file of the fourier-transformed data"
    write(*,*) "        if not referred, 'signal.dat' ist used"
    write(*,*) "        overrides a previous --no-signal option"
    write(*,*) "--no-signal"
    write(*,*) "        the fourier-transformed data is not given out"
    write(*,*) "        overrides a previous -s of --signal option"
    write(*,*) "-K <integer>"
    write(*,*) "        the number of points which are used to calculate the energies and amplitudes"
    write(*,*) "        it should be an integer between 50 and 200, but it mustn't be larger than the number of input data points"
    write(*,*) "        if not referred, 100 is used"
    write(*,*) "-g <float>, --genau <float>"
    write(*,*) "        the precision which is used to identify exceptional points"
    write(*,*) "        it should be an floating point number between 1e-3 and 1e-7,"
    write(*,*) "            but it mustn't be greater than zero"
    write(*,*) "        if not referred, 1e-5 is used"

  end subroutine print_help

end program hi
