! gfortran -o normal_state_E.out normal_state_E.f90 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
program main
    implicit none
    !parameter
    integer, parameter :: N = 400, allhop = 800, loop_num = 300
    double precision, parameter :: pi = 4*atan(1.), hsq = 0., max_V = 3.0 ! move V
    character(7), parameter :: hop_file = "hop.txt"
    character(11), parameter :: E_file = "eigvals.txt"
    !variable
    integer :: unit_write_result, loop, index
    double precision :: V = 0., h_z = sqrt(hsq)
    complex(kind(0d0)) :: Hamiltonian(2 * N, 2 * N)
    !Density of state of zero energy
    double precision :: Es(2 * N, loop_num)
    ! for ZHEEVD of lapack
    integer :: INFO, IWORK(1)
    double precision :: W(2 * N), RWORK(2 * N)
    complex(kind(0d0)) :: WORK(2 * N + 1)

    do loop = 1, loop_num
        V = max_V * (loop - 1) / (loop_num - 1)
        call make_Hamiltonian()
        call ZHEEVD('N', 'U', 2 * N, Hamiltonian, 2 * N, W, WORK, 2 * N + 1,&
        RWORK, 2 * N, IWORK, 1, INFO)

        do index = 1, 2 * N
            Es(index, loop) = W(index)
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

    subroutine write_files()
        integer :: i

        open(newunit=unit_write_result, file=E_file)
            do i = 1, 2 * N
                write(unit_write_result, *) Es(i, 1:loop_num)
            end do
        close(unit_write_result)
    end subroutine write_files
end program main