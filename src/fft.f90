module fft
  use myKinds
  use fftw3constants
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! this module does the FFT using the fftw library
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! make the important functions public
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! subroutine to ensure right array lengths and call the FFT-routine
  public :: fft_fourier

contains

  ! subroutine to ensure right array lengths and call the FFT-routine
  subroutine fft_fourier(min,max,N,K,f,fourier,tau)
    implicit none

    real(dp),    intent(in)  :: min          ! minimum of energy window
    real(dp),    intent(in)  :: max          ! maximum of energy window
    integer,     intent(in)  :: N            ! length of input arrays
    ! length of signal to be transformed
    integer,     intent(in)  :: K
    ! input signal
    complex(dp), dimension(N), intent(in)  :: f
    ! output signal (FFT-transformed)
    complex(dp), dimension(N), intent(out) :: fourier
    ! time step tau
    real(dp),    intent(out) :: tau

    ! signals in the right length
    complex(dp), dimension(K) :: fK, fourierK

    if(K > N) then
       write(0,*) "error in fft: length of array is smaller than&
            & number of data points"
       call exit(2)
    end if

    fK(1:K) = f(1:K)

    call fft_fourierk(min,max,K,fK,fourierK,tau)

    fourier(1:K) = fourierK(1:K)
    if(N > K) fourier(K+1:N) = (0d0,0d0)

  end subroutine fft_fourier

  ! subroutine to do the FFT
  subroutine fft_fourierk(min,max,K,f,fourier,tau)
    implicit none

    real(dp),    intent(in)  :: min, max
    integer,     intent(in)  :: K
    complex(dp), dimension(K), intent(in)  :: f
    complex(dp), dimension(K), intent(out) :: fourier
    real(dp),    intent(out) :: tau

    complex(dp), dimension(K) :: f_work

    integer(dp) :: plan

    tau = 8d0*atan(1d0)/(max-min)

    ! fft
    call dfftw_plan_dft_1d(plan, K, f, f_work, FFTW_FORWARD,&
         & FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)

    ! ensure right prefactor
    fourier(1:K) = (max-min)/(K*8d0*atan(1d0)) * f_work(1:K)
    
  end subroutine fft_fourierk
  
end module fft
