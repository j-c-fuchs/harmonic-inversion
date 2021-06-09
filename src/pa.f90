module pa
  use myKinds
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! this module contains routines to solve the equations via Pade-
!!! approximation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! some general comments:
!!!   - for degree of degeneration less than or equal to 3 the explicit
!!!     formulas from the paper are used, see Eqn. (26) to (28))
!!!   - the implementation uses some lapack routines to solve standard
!!!     algebraic problems
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! make the important functions public
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! calculates the resonaces (the results are not transformed, i.e.
  ! they are in the form of d^ in the paper)
  public :: pa_all
  ! transform results to normal form (i.e. d in the paper)
  public :: pa_skal
  ! transform results to fourier-like form (i.e. d~ in the paper, see
  ! Eqn. (36))
  public :: pa_skal_fourier

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! public functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! calculates the resonaces (not transformed)
  subroutine pa_all(N, K, tau, prec_deg, MEAN, GEOMETRIC,&
       & c, z, z_dod, w, d)
    implicit none

    ! for the length of input arrays
    integer,  intent(in) :: N
    ! parameter K (number of resonances to extract)
    integer,  intent(in) :: K
    ! parameter tau
    real(dp), intent(in) :: tau
    ! maximum distance between to resonances to treat them as
    ! degenerated  (parameter delta in the paper, p. 7)
    real(dp), intent(in), optional :: prec_deg

    ! compute mean of the frequencies?
    ! (compute means of the z's of degenerated resonances)
    logical,  intent(in) :: MEAN
    ! compute geometric mean instead of the arithmetic of the
    ! frequencies?
    ! (taking the geometric mean means taking the arithmetic mean of
    ! the corresponding w's)
    logical,  intent(in) :: GEOMETRIC

    ! signal c to analyse (twice as long as other arrays!)
    complex(dp), dimension(2*N), intent(in)  :: c

    ! wanted z's, w's and d's
    complex(dp), dimension(N), intent(out) :: z,w,d
    ! corresponding degree of degeneration
    integer, dimension(N), intent(out) :: z_dod
    
    ! coefficients of the polynomials P and Q
    complex(dp), dimension(K) :: a,b

    ! calculate coefficients of the first polynomial (Q)
    call pa_coeffa(N, K, c, a)
    ! calculate coefficients of the second polynomial (P)
    call pa_coeffb(N, K, a, c, b)
    ! calculate the roots of the polynomial P (i.e. the wanted z's) via
    ! the diagonalisation of the Hessenberg Matrix from Eqn. (10)
    call pa_paramz(N, K, a, z)
    ! compare the calculated roots to find the (possibly) degenerated
    ! resonances
    call pa_compare(N, K, prec_deg, MEAN, GEOMETRIC, z, z_dod)
    ! calculate energies from frequencies
    call pa_energies(N, K, tau, z, w)
    ! calculate amplitudes
    call pa_amplit(N, K, a, b, z, z_dod, d)

  end subroutine pa_all

  ! transform results to normal form
  subroutine pa_skal(N, K, min, max, tau, z_dod, w, d)
    implicit none

    integer,  intent(in)  :: N,K
    real(dp), intent(in)  :: min,max,tau
    integer,     dimension(N), intent(in)    :: z_dod
    complex(dp), dimension(N), intent(inout) :: w,d

    integer :: i,j

    ! scale amplitudes
    i = 1
    do while(i <= K)
       ! degree of degeneration is 1
       if(z_dod(i) == 1) then
          d(i) = d(i) * (0_dp,1_dp)
          i = i + 1
          cycle
       end if
       ! degree of degeneration is 2
       if(z_dod(i) == 2) then
          d(i) = d(i) * (-1_dp/tau)
          d(i+1) = d(i+1) * (0_dp,1_dp)
          i = i + 2
          cycle
       end if
       ! degree of degeneration is 3
       if(z_dod(i) == 3) then
          d(i) = d(i) * (0_dp,-2_dp)/(tau*tau)
          d(i+1) = d(i+1) * (-1_dp/tau)
          d(i+2) = d(i+2) * (0_dp,1_dp)
          i = i + 3
          cycle
       end if
       ! degree of degeneration is higher than 3
       if(z_dod(i) > 3) then
          do j = 0, z_dod(i)-1
             d(i+j) = (0_dp,1_dp) * d(i+j) * pa_math_faculty(j)&
                  & / (((0_dp,-1_dp)*tau)**j)
          end do
          i = i + z_dod(i)
          cycle
       end if
       ! avoid endless loops if z_dod < 1
       i = i + 1
    end do

    ! shift frequencies
    do i=1, K
       if(dble(w(i)) < 0) w(i) = w(i) + (max-min)
       w(i) = w(i) + min
    end do

  end subroutine pa_skal

  ! transform results to fourier-like form
  subroutine pa_skal_fourier(N, K, min, max, tau, z_dod, w, d)
    implicit none

    integer,  intent(in)  :: N,K
    real(dp), intent(in)  :: min,max,tau
    integer,     dimension(N), intent(in)    :: z_dod
    complex(dp), dimension(N), intent(inout) :: w,d

    integer :: i,j

    ! scale amplitudes
    i = 1
    do while(i <= K)
       ! degree of degeneration is 1
       if(z_dod(i) == 1) then
          i = i + 1
          cycle
       end if
       ! degree of degeneration is 2
       if(z_dod(i) == 2) then
          d(i) = d(i) / tau
          i = i + 2
          cycle
       end if
       ! degree of degeneration is 3
       if(z_dod(i) == 3) then
          d(i)   = d(i)   / (tau*tau)
          d(i+1) = d(i+1) / tau
          i = i + 3
          cycle
       end if
       ! degree of degeneration is higher than 3
       if(z_dod(i) > 3) then
          do j = i, i+z_dod(i)-1
             d(j) = d(j) / tau**(i+z_dod(i)-j)
          end do
          i = i + z_dod(i)
          cycle
       end if
       ! avoid endless loops if z_dod < 1
       i = i + 1
    end do

    ! frequencies needn't to be shifted in this case

  end subroutine pa_skal_fourier

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! private functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! calculate coefficients of the first polynomial
  subroutine pa_coeffa(N,K,c,a)
    implicit none

    ! This subroutine uses the routine zgesv from the lapack-library to
    ! solve the system of linear equations which is built for general
    ! systems of linear equations.
    ! You might also have a look at the routines
    !   - zgerfs (does iterative refinement to improve the result)
    !   - zgesvx (zgesv with zgerfs built in)
    !   - zsysv  (routine for complex symmetric matrices)
    !   - zsysvx (zsysv with zgerfs built in)
    !   - zgelss (routine for underdetermined systems) or
    !   - zgelsy (routine for underdetermined, complex symmetric
    !             systems).

    integer, intent(in)  :: N,K
    complex(dp), dimension(2*N), intent(in)  :: c
    complex(dp), dimension(K), intent(out) :: a

    ! matrix for the system of linear equations
    complex(dp), dimension(K,K) :: M

    integer :: i,j,INFO

    ! for the lapack routine which solves the system of linear equation
    integer, dimension(K) :: IPIV

    ! prepare the matrix and right hand side
    forall(j = 1:K, i = 1:K) M(i,j) = c(i+j)
    a(1:K) = c(1:K)

    ! solve the system of linear equations
    call zgesv(K, 1, M, K, IPIV, a, K, INFO)
    ! treat errors
    if (INFO /= 0) then
       write(0,*) "pa_coeffa: error",INFO,"in function 'zgesv'&
            & to compute the Polynom P"
       call exit(2)
    end if

  end subroutine pa_coeffa

  ! calculate coefficients of the second polynomial (P)
  subroutine pa_coeffb(N,K,a,c,b)
    implicit none

    integer, intent(in)  :: N,K
    complex(dp), dimension(K), intent(in)   :: a
    complex(dp), dimension(2*N), intent(in) :: c
    complex(dp), dimension(K), intent(out)  :: b

    integer :: i,j

    do i=1, K
       b(i) = (0_dp,0_dp)
       do j=0, K-i
          b(i) = b(i) + a(i+j)*c(j+1)
       end do
    end do

  end subroutine pa_coeffb

  ! calculate the roots of the polynomial P (i.e. the wanted z's)
  subroutine pa_paramz(N,K,a,z)
    implicit none

    ! This subroutine uses the routine zhseqr from the lapack-library
    ! to compute the eigenvalues of the Hessenberg matrix (see Eqn.
    ! (10) in the paper).
    ! Alternatives could be the routines
    !   - zgees (for general matrices) or
    !   - zgeev (ditto).

    integer, intent(in)  :: N,K
    complex(dp), dimension(K), intent(in)  :: a
    complex(dp), dimension(N), intent(out) :: z

    ! matrix for the calculation of the eigenvalues
    complex(dp), dimension(K,K) :: M

    integer :: i

    complex(dp), dimension(K) :: eigenvalues

    ! for the lapack routine which computes the eigenvalues
    integer :: INFO
    integer :: ILO
    integer :: IHI
    complex(dp), dimension(K,K) :: SCHUR
    complex(dp), dimension(2*K*K) :: WORK
    integer :: LWORK

    ILO = 1
    IHI = K
    LWORK = 2*K*K

    ! prepare the matrix
    M(1:K,1:K) = (0_dp,0_dp)
    do i=2, K
       M(i,i-1) = (1_dp,0_dp)
    end do
    do i=1, K-1
       M(1,i) = -a(K-i)/a(K)
    end do
    M(1,K) = (1_dp,0_dp)/a(K)

    ILO = 1
    IHI = K

    ! compute the eigenvalues
    call zhseqr('E', 'N', K, ILO, IHI, M, K, eigenvalues, schur, K,&
         & WORK, LWORK, INFO)
    ! treat errors
    if(INFO /= 0) then
       write(0,*) "error", INFO, "in function 'zhseqr'&
            & to compute the eigenvalues resp. zeros"
       call exit(2)
    end if

    z(1:K) = eigenvalues(1:K)
    z(K+1:N) = (0_dp,0_dp)

  end subroutine pa_paramz

  ! compare the calculated roots to find the (possibly) degenerated
  ! resonances
  subroutine pa_compare(N, K, prec_deg, MEAN, GEOMETRIC, z, z_dod)
    implicit none

    integer, intent(in) :: N,K
    real(dp), intent(in) :: prec_deg
    logical, intent(in) :: MEAN
    logical, intent(in) :: GEOMETRIC

    complex(dp), dimension(N), intent(inout) :: z
    integer, dimension(N), intent(out) :: z_dod
    
    integer :: i

    complex(dp), dimension(K) :: a
    integer, dimension(K) :: v

    a(1:K) = z(1:K)

    ! pa_compareL compares the values
    call pa_compareL(K, prec_deg, MEAN, GEOMETRIC, a, v)

    z(1:K) = a(1:K)
    z_dod(1:K) = v(1:K)

    ! avoid endless loops
    do i = K+1, N
       z_dod(i) = N-i
    end do

  end subroutine pa_compare

  ! pa_compareL compares a value with all following values
  ! NOTE: speed can be further optimized
  subroutine pa_compareL(K, prec_deg, MEAN, GEOMETRIC, a, v)
    implicit none

    integer, intent(in) :: K
    real(dp), intent(in) :: prec_deg
    logical, intent(in) :: MEAN
    logical, intent(in) :: GEOMETRIC
    complex(dp), dimension(K), intent(inout) :: a
    integer,     dimension(K), intent(out)   :: v

    integer :: i,j
    integer :: dod
    complex(dp) :: mittel

    complex(dp), dimension(K) :: swap
    integer,     dimension(K) :: number
    integer,     dimension(K) :: number_dod
    integer :: NO_count
    integer :: EV_count

    ! compare values and get degree of degeneration
    number(1:K) = 0
    NO_count = 0
    do i=1, K
       if(number(i) > 0) cycle
       NO_count = NO_count + 1
       number(i) = NO_count
       dod = 1
       do j=i+1, K
          if(number(j) > 0) cycle
          if(abs(a(i)-a(j)) > prec_deg) cycle
          dod = dod + 1
          number(j) = NO_count     ! mark (possibly) degenerated values
       end do
       number_dod(NO_count) = dod
    end do
     
    ! sort values via array swap
    EV_count = 0
    do i=1, NO_count
       do j=1, K
          if(number(j) /= i) cycle
          EV_count = EV_count + 1
          swap(EV_count) = a(j)
          v(EV_count) = number_dod(i)
       end do
    end do

    a(1:K) = swap(1:K)

    if(MEAN) then
       ! compute the mean of degenerated values
       if(GEOMETRIC) then
          ! use the geometric mean
          i = 1
          do while(i <= K)
             if(v(i) < 2) then
                i = i + 1
                cycle
             end if
             mittel = (1_dp,0_dp)
             do j = i, i+v(i)-1
                mittel = mittel * a(j)
             end do
             mittel = exp( log(mittel)/v(i) )
             do j = i, i+v(i)-1
                a(j) = mittel
             end do
             i = i + v(i)
          end do
       else
          ! use the arithmetic mean
          i = 1
          do while(i <= K)
             if(v(i) < 2) then
                i = i + 1
                cycle
             end if
             mittel = (0_dp,0_dp)
             do j = i, i+v(i)-1
                mittel = mittel + a(j)
             end do
             mittel = mittel / v(i)
             do j = i, i+v(i)-1
                a(j) = mittel
             end do
             i = i + v(i)
          end do
       end if
    else
       ! don't compute the mean
       ! instead set the following values to the first one
       ! (so that all have the same value)
       i = 1
       do while(i <= K)
          if(v(i) < 2) then
             i = i + 1
             cycle
          end if
          do j = i+1, i+v(i)-1
             a(j) = a(i)
          end do
          i = i + v(i)
       end do
    end if
    
    if(EV_count /= K) then
       write(0,*) "error while sorting the eigenvalues:&
            & number of eigenvalues (",EV_count,") do not match", K
       call exit(2)
    end if
    
  end subroutine pa_compareL

  ! calculate energies from frequencies
  subroutine pa_energies(N,K,tau,z,w)
    implicit none

    integer, intent(in) :: N,K
    real(dp),    intent(in) :: tau
    complex(dp), dimension(N), intent(in)  :: z
    complex(dp), dimension(N), intent(out) :: w

    complex(dp) :: tauI

    tauI = (0_dp,1_dp) / tau
    w(1:K) = tauI * log(z(1:K))

  end subroutine pa_energies

  ! calculate amplitudes
  subroutine pa_amplit(N,K,a,b,z,z_dod,amplit)
    implicit none

    integer, intent(in)  :: N, K
    complex(dp), dimension(K), intent(in)  :: a, b
    complex(dp), dimension(N), intent(in)  :: z
    integer, dimension(N), intent(in)  :: z_dod
    complex(dp), dimension(N), intent(out) :: amplit

    integer :: i,j,ji,jj
    complex(dp) :: sum

    complex(dp) :: Q1 = (0_dp,0_dp)
    complex(dp) :: Q2 = (0_dp,0_dp)
    complex(dp) :: Q3 = (0_dp,0_dp)
    complex(dp) :: Q4 = (0_dp,0_dp)

    i = 1
    do while(i <= K)

       ! do not divide by zero
       if(z(i) == (0_dp,0_dp)) then
          write(0,*) "pa_amplit: eigenvalue",i,"is zero"
          amplit(i) = (0_dp,0_dp)
          i = i + 1
          cycle
       end if

       ! degree of degeneration is 1
       if(z_dod(i) == 1) then
          Q1 = pa_Q1(K,a,z(i))
          ! do not divide by zero
          if(Q1 == (0_dp,0_dp)) then
             write(0,*) "error: point number",i,"may be an EP&
                  & (because Q' = 0) but wasn't identified as this;&
                  & no amplitude computed"
             amplit(i) = (0_dp,0_dp)
             i = i + 1
             cycle
          end if
          amplit(i) = pa_P0(K,b,z(i)) / (z(i)*Q1)
          i = i + 1
          cycle
       end if

       ! degree of degeneration is 2
       if(z_dod(i) == 2) then
          Q2 = pa_Q2(K,a,z(i))
          ! do not divide by zero
          if(Q2 == (0_dp,0_dp)) then
             write(0,*) "error: point number",i,"may be an EP of higher&
                  & order than 2 (because Q'' = 0) but wasn't&
                  & identified as this; no amplitude computed"
             amplit(i) = (0_dp,0_dp)
             amplit(i+1) = (0_dp,0_dp)
             i = i + 2
             cycle
          end if
          amplit(i) = pa_P0(K,b,z(i)) / ( z(i)*z(i) * Q2 )
          sum = z(i)*pa_Q3(K,a,z(i)) / Q2 + 1_dp
          amplit(i+1) = pa_P1(K,b,z(i)) / (z(i)*Q2) - amplit(i)*sum
          i = i + 2
          cycle
       end if

       ! degree of degeneration is 3
       if(z_dod(i) == 3) then
          Q3 = pa_Q3(K,a,z(i))
          Q4 = pa_Q4(K,a,z(i))
          ! do not divide by zero
          if(Q3 == (0_dp,0_dp)) then
             write(0,*) "error: point number",i,"may be an EP of higher&
                  & order than 3 (because Q'' = 0) but wasn't&
                  & identified as this; no amplitude computed"
             amplit(i) = (0_dp,0_dp)
             amplit(i+1) = (0_dp,0_dp)
             amplit(i+2) = (0_dp,0_dp)
             i = i + 3
             cycle
          end if
          amplit(i) = pa_P0(K,b,z(i)) / (2_dp*z(i)*z(i)*z(i)*Q3)
          sum = 2_dp*z(i)*Q4 / Q3 + 3_dp
          amplit(i+1) = pa_P1(K,b,z(i)) / (z(i)*z(i)*Q3) - amplit(i)*sum
          sum = z(i)*Q4 / Q3 + 1_dp
          amplit(i+2) = pa_P2(K,b,z(i)) / (z(i)*Q3) - amplit(i+1)*sum
          sum = 2_dp*z(i)*z(i)*pa_Q5(K,a,z(i)) / Q3 &
               &+ 3_dp*z(i)*Q4 / Q3 + 1_dp
          amplit(i+2) = amplit(i+2) - amplit(i)*sum
          i = i + 3
          cycle
       end if

       ! degree of degeneration is higher than 9
       ! do nothing because of long computational time
       if(z_dod(i) > 9) then
          write(*,*) "warning: point number",i,"is a EP of higher order&
               & than 9; no amplitudes computed because of long&
               & computation time"
          do j = z_dod(i), 1, -1
             amplit(i+z_dod(i)-j) = (0_dp,0_dp)
          end do
          i = i + z_dod(i)
          cycle
       end if

       ! degree of degeneration is higher than 3
       if(z_dod(i) > 3) then
          ! do not divide by zero
          if(pa_Q(K,a,z(i),z_dod(i)) == (0_dp,0_dp)) then
             write(0,*) "error: point number",i,"may be an EP of higher&
                  & order than",K," (because the ",K,"-th derivative&
                  & of Q is zero) but wasn't identified as this;&
                  & no amplitude computed"
             do j = i, i+z_dod(i)-1
                amplit(j) = (0_dp,0_dp)
             end do
             i = i + z_dod(i)
             cycle
          end if
          do j = z_dod(i), 1, -1
             amplit(i+z_dod(i)-j) = pa_P(K,b,z(i),z_dod(i)-j)&
                  & / (z(i)**j)
             do ji = z_dod(i), j+1, -1
                sum = (0_dp,0_dp)
                do jj = 0, ji-j
                   sum = sum + pa_math_stirling2_part(ji,j+jj)&
                        & * (z(i)**jj) * pa_Q(K,a,z(i),z_dod(i)+jj)&
                        & / (j+jj)
                end do
                amplit(i+z_dod(i)-j) = amplit(i+z_dod(i)-j)&
                     & - amplit(i+z_dod(i)-ji) * sum
             end do
             amplit(i+z_dod(i)-j) = amplit(i+z_dod(i)-j)&
                  & / ( pa_math_faculty(j-1) * pa_Q(K,a,z(i),z_dod(i)) )
          end do
          i = i + z_dod(i)
          cycle
       end if
       
       ! avoid endless loops
       i = i + 1
       
    end do

  end subroutine pa_amplit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! help functions (polynomials)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! function to compute the d-th derivative of the polynomial P at x
  ! devided by n factorial
  complex(dp) function pa_P(K, b, x, d)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: b         ! coefficients b
    complex(dp), intent(in) :: x   ! point
    integer, intent(in) :: d       ! order of derivation

    complex(dp), dimension(0:K) :: array

    array(0) = (0_dp,0_dp)
    array(1:K) = b(1:K)
    
    pa_P = pa_math_polynom_taylor_cdp(K, array, x, d)

  end function pa_P

  ! function to compute the d-th derivative of the polynomial Q at x
  ! devided by n factorial
  complex(dp) function pa_Q(K, a, x, d)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: a         ! coefficients a
    complex(dp), intent(in) :: x   ! point
    integer, intent(in) :: d       ! order of derivation

    complex(dp), dimension(0:K) :: array

    array(0) = (-1_dp,0_dp)
    array(1:K) = a(1:K)
    
    pa_Q = pa_math_polynom_taylor_cdp(K, array, x, d)

  end function pa_Q

  ! some of these derivatives are implemented directly
  ! (pa_P0(K, b, x) = pa_P(K, b, x, 0) and so on)

  ! function to compute the value of the polynomial P at x
  complex(dp) function pa_P0(K, b, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: b         ! coefficients b
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_P0 = b(K)
    do i = 1, K-1
       pa_P0 = pa_P0 * x + b(K-i)
    end do
    pa_P0 = pa_P0 * x

  end function pa_P0

  ! function to compute the first derivative of the polynomial P at x
  complex(dp) function pa_P1(K, b, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: b         ! coefficients b
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_P1 = real(K, dp) * b(K)
    do i = 1, K-1
       pa_P1 = pa_P1 * x + real(K-i, dp) * b(K-i)
    end do

  end function pa_P1

  ! function to compute the second derivative of the polynomial P at x
  ! devided by 2
  complex(dp) function pa_P2(K, b, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: b         ! coefficients b
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_P2 = real(K, dp) * real(K-1, dp) * b(K)
    do i = 1, K-2
       pa_P2 = pa_P2 * x + real(K-i, dp) * real(K-i-1, dp) * b(K-i)
    end do
    pa_P2 = pa_P2 / 2_dp

  end function pa_P2

  ! function to compute the value of the polynomial Q at x
  complex(dp) function pa_Q0(K, a, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: a         ! coefficients a
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_Q0 = a(K)
    do i = 1, K-1
       pa_Q0 = pa_Q0 * x + a(K-i)
    end do
    pa_Q0 = pa_Q0 * x - 1_dp

  end function pa_Q0

  ! function to compute the first derivative of the polynomial Q at x
  complex(dp) function pa_Q1(K, a, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: a         ! coefficients a
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_Q1 = real(K, dp) * a(K)
    do i = 1, K-1
       pa_Q1 = pa_Q1 * x + real(K-i, dp) * a(K-i)
    end do

  end function pa_Q1

  ! function to compute the second derivative of the polynomial Q at x
  ! devided by 2
  complex(dp) function pa_Q2(K, a, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: a         ! coefficients a
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_Q2 = real(K, dp) * real(K-1, dp) * a(K)
    do i = 1, K-2
       pa_Q2 = pa_Q2 * x + real(K-i, dp) * real(K-i-1, dp) * a(K-i)
    end do
    pa_Q2 = pa_Q2 / 2_dp

  end function pa_Q2

  ! function to compute the third derivative of the polynomial Q at x
  ! devided by 6
  complex(dp) function pa_Q3(K, a, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: a         ! coefficients a
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_Q3 = real(K, dp) * real(K-1, dp) * real(K-2, dp) * a(K)
    do i = 1, K-3
       pa_Q3 = pa_Q3 * x + real(K-i, dp) * real(K-i-1, dp)&
            & * real(K-i-2, dp) * a(K-i)
    end do
    pa_Q3 = pa_Q3 / 6_dp

  end function pa_Q3

  ! function to compute the fourth derivative of the polynomial Q at x
  ! devided by 24
  complex(dp) function pa_Q4(K, a, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: a         ! coefficients a
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_Q4 = real(K, dp) * real(K-1, dp) * real(K-2, dp) * real(K-3, dp)&
         & * a(K)
    do i = 1, K-4
       pa_Q4 = pa_Q4 * x + real(K-i, dp) * real(K-i-1, dp)&
            & * real(K-i-2, dp) * real(K-i-3, dp) * a(K-i)
    end do
    pa_Q4 = pa_Q4 / 24_dp

  end function pa_Q4

  ! function to compute the fifth derivative of the polynomial Q at x
  ! devided by 120
  complex(dp) function pa_Q5(K, a, x)
    implicit none

    integer, intent(in) :: K       ! number of coefficients
    complex(dp), dimension(K), intent(in) :: a         ! coefficients a
    complex(dp), intent(in) :: x   ! point

    integer :: i

    pa_Q5 = real(K, dp) * real(K-1, dp) * real(K-2, dp) * real(K-3, dp)&
         & * real(K-4, dp) * a(K)
    do i = 1, K-5
       pa_Q5 = pa_Q5 * x + real(K-i, dp) * real(K-i-1, dp)&
            & * real(K-i-2, dp) * real(K-i-3, dp) * real(K-i-4, dp)&
            & * a(K-i)
    end do
    pa_Q5 = pa_Q5 / 120_dp

  end function pa_Q5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! mathematical help functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! faculty
  real(dp) function pa_math_faculty(n)
    implicit none

    integer, intent(in) :: n

    integer :: i

    pa_math_faculty = 1_dp

    if(n < 0) then
       write(0,*) "pa_math_faculty: input variable < 0"
       call exit(1)
    elseif(n == 0 .or. n == 1) then
       return
    elseif(n > 32) then
       pa_math_faculty = dgamma(n+1d0)
    else
       do i = 2, n
          pa_math_faculty = pa_math_faculty * i
       end do
    end if

  end function pa_math_faculty

  ! stirling numbers of second kind
  real(dp) function pa_math_stirling2_part(n,k)
    implicit none

    integer, intent(in) :: n,k

    integer :: i
    real(dp) :: out
    
    pa_math_stirling2_part = 0_dp

    if(k > n .or. k < 0) return

    out = 0_dp
    do i = 0, k
       if(mod(i,2) == 0) then
          out = out + pa_math_binom(k,i)*(real(i, dp)**n)
       else
          out = out - pa_math_binom(k,i)*(real(i, dp)**n)
       end if
    end do
    pa_math_stirling2_part = abs(out)
    
  end function pa_math_stirling2_part

  ! binomial coeffizients
  real(dp) function pa_math_binom(n,k)
    implicit none

    integer, intent(in) :: n,k

    integer :: i,j,l

    if(k < 0 .or. k > n) then
       pa_math_binom = 0_dp
    elseif(k == 0 .or. k == n) then
       pa_math_binom = 1_dp
    elseif(k > 32) then
       pa_math_binom = dgamma(n+1d0)/(dgamma(k+1d0)+dgamma(n-k+1d0))
    else
       l = k
       if(2*k > n) l = n-k
       j = n-l+1
       do i=2,l
          j = j * (n-l+i)
          j = j / i
       end do
       pa_math_binom = real(j, dp)
    end if

  end function pa_math_binom

  ! d-th derivative of polynomial with coeffizients a0,a1,...,aN at x
  real(dp) function pa_math_polynom_taylor(N,a,x,d)
    implicit none

    integer, intent(in) :: N
    real(dp), dimension(0:N), intent(in) :: a
    real(dp), intent(in) :: x
    integer, intent(in) :: d

    integer :: i,j
    real(dp), dimension(0:N,0:N+1) :: b

    if(d < 0) then
       write(0,*) "mymath_polynom: order of derivative must be&
            & non-negative"
       call exit(1)
    end if

    b(0:N,0) = a(0:N)
    b(N,1:N+1) = a(N)
    do i = 1, d+1
       do j = N-1, i-1, -1
          b(j,i) = b(j,i-1) + b(j+1,i)*x
       end do
    end do
    pa_math_polynom_taylor = b(d,d+1)
    
  end function pa_math_polynom_taylor

  ! d-th derivative of polynomial with coeffizients a0,a1,...,aN at x
  ! (complex version)
  complex(dp) function pa_math_polynom_taylor_cdp(N,a,x,d)
    implicit none

    integer, intent(in) :: N
    complex(dp), dimension(0:N), intent(in) :: a
    complex(dp), intent(in) :: x
    integer, intent(in) :: d

    integer :: i,j
    complex(dp), dimension(0:N,0:N+1) :: b

    if(d < 0) then
       write(0,*) "pa_math_polynom_cdp: order of derivative must be&
            & non-negative"
       call exit(1)
    end if

    b(0:N,0) = a(0:N)
    b(N,1:N+1) = a(N)
    do i = 1, d+1
       do j = N-1, i-1, -1
          b(j,i) = b(j,i-1) + b(j+1,i)*x
       end do
    end do
    pa_math_polynom_taylor_cdp = b(d,d+1)
    
  end function pa_math_polynom_taylor_cdp

end module pa
