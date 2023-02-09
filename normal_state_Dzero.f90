! gfortran -o normal_state_Dzero.out normal_state_Dzero.f90 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
program main
    implicit none
    !parameter
    integer, parameter :: N = 400, allhop = 800, points_num = 10**4 + 1, loop_num1 = 100, loop_num2 = 300
    double precision, parameter :: pi = 4*atan(1.), gamma = 0.05, max_hsq = 1.0, max_V = 3.0
    double precision, parameter :: width = 10.0 ! エネルギー固有値の幅(-width < E < width)
    character(7), parameter :: hop_file = "hop.txt"
    character(12), parameter :: D0_file = "DOS_zero.txt"
    !variable
    integer :: unit_write_result, loop1, loop2
    double precision :: hsq = 0., V = 0., h_z
    complex(kind(0d0)) :: Hamiltonian(2 * N, 2 * N)
    !Density of state of zero energy
    double precision :: D_0(loop_num1, loop_num2) = 0
    ! for ZHEEVD of lapack
    integer :: INFO, IWORK(1)
    double precision :: W(2 * N), RWORK(2 * N)
    complex(kind(0d0)) :: WORK(2 * N + 1)

    do loop1 = 1, loop_num1
        hsq = max_hsq * (loop1 - 1) / (loop_num1 - 1)
        h_z = sqrt(hsq)
        do loop2 = 1, loop_num2
            V = max_V * (loop2 - 1) / (loop_num2 - 1)
            call make_Hamiltonian()
            call ZHEEVD('N', 'U', 2 * N, Hamiltonian, 2 * N, W, WORK, 2 * N + 1,&
            RWORK, 2 * N, IWORK, 1, INFO)

            D_0(loop1, loop2) = DOS_zero()
        end do
    end do

    call write_files()
contains
    subroutine make_Hamiltonian()
        integer :: i, j, k, l, m
        double precision :: x, y, r
        complex(kind(0d0)) :: H_ij

        Hamiltonian = 0

        ! diagonal elements
        do k = 1, N
            ! spin up
            i = 2 * k - 1
            j = i
            H_ij = h_z

            Hamiltonian(i, j) = H_ij

            ! spin down
            i = 2 * k
            j = i
            H_ij = -h_z !

            Hamiltonian(i, j) = H_ij
        end do

        ! hopping elements
        open(newunit = unit_write_result, file = hop_file)
        do k = 1, allhop
            read(unit_write_result, *) l, m
            l = l + 1
            m = m + 1

            ! spin up
            i = 2 * l - 1
            j = 2 * m - 1
            H_ij = -1.

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = H_ij

            ! spin down
            i = 2 * l
            j = 2 * m
            H_ij = -1.

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = H_ij
        end do
        close (unit_write_result)

        !spin-orbit coupling
        open (newunit = unit_write_result, file=hop_file)
        do k = 1, allhop
            read (unit_write_result, *) l, m, x, y
            l = l + 1
            m = m + 1

            r = sqrt(x**2 + y**2)
            x = x/r
            y = y/r

            !iup, jdown
            i = 2*l - 1
            j = 2*m
            H_ij = V*cmplx(-x, y, kind(0d0))

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = conjg(H_ij)

            !idown, jup
            i = 2*l
            j = 2*m - 1
            H_ij = V*cmplx(x, y, kind(0d0))

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = conjg(H_ij)
        end do
        close(unit_write_result)
    end subroutine make_Hamiltonian

    function Lorentzian(x, eigen_energy) result(result)
        double precision :: eigen_energy
        double precision :: x(points_num), result(points_num)

        result = gamma / (pi * (gamma**2 + (x - eigen_energy)**2))
    end function

    function DOS_zero() result(result)
        integer :: i
        double precision, parameter :: x(points_num) = [(-width + 2 * width * (i - 1) / (points_num - 1), i = 1, points_num)]
        double precision :: DOS(points_num)
        double precision :: result

        DOS = 0

        do i = 1, 2 * N
            DOS = DOS + Lorentzian(x, W(i))
        end do

        result = DOS(points_num/2 + 1)
    end function DOS_zero

    subroutine write_files()
        integer :: i

        open(newunit=unit_write_result, file=D0_file)
            do i = 1, loop_num1
                write(unit_write_result, *) D_0(i, 1:loop_num2)
            end do
        close(unit_write_result)
    end subroutine write_files
end program main
