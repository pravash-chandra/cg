program main
        use mpi
        implicit none

        integer, dimension(7) :: n_values
        integer :: n, test_idx
        real, allocatable :: A(:,:), b(:), x(:)
        integer :: i, j
        double precision :: start_time, end_time, elapsed_time

        ! MPI variables
        integer :: rank, nprocs, ierr
        integer :: rows_local, row_start, row_end
        real, allocatable :: A_local(:,:), Ap_local(:)

        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

        n_values = [5000, 10000, 15000, 20000, 25000, 30000, 35000]

        if (rank == 0) then
                print *, "====== CG Solver Benchmarking (MPI Parallel) ======"
                print *, "Processes:", nprocs
                print *, ""
                print '(A15, A15)', "Matrix Size", "Time (sec)"
                print *, "-------------------------------------------"
        end if

        do test_idx = 1, 7
                n = n_values(test_idx)

                ! All ranks allocate full A and b (simple approach)
                if (allocated(A)) deallocate(A)
                if (allocated(b)) deallocate(b)
                if (allocated(x)) deallocate(x)
                allocate(A(n,n), b(n), x(n))

                ! Build matrix and rhs on all ranks
                do i = 1, n
                        do j = 1, n
                                A(i,j) = 1.0 / (abs(i-j) + 1.0)
                        end do
                end do
                do i = 1, n
                        A(i,i) = A(i,i) + real(n)
                end do
                b = 1.0

                ! Each rank owns a contiguous block of rows
                rows_local = n / nprocs
                row_start  = rank * rows_local + 1
                row_end    = row_start + rows_local - 1
                if (rank == nprocs - 1) row_end = n   ! last rank takes remainder

                allocate(A_local(row_end - row_start + 1, n))
                allocate(Ap_local(row_end - row_start + 1))
                A_local = A(row_start:row_end, :)

                start_time = MPI_Wtime()
                x = cg_mpi(A_local, b, n, row_start, row_end, rank, nprocs)
                end_time   = MPI_Wtime()
                elapsed_time = end_time - start_time

                if (rank == 0) then
                        print '(I15, F15.6)', n, elapsed_time
                end if

                deallocate(A_local, Ap_local)
        end do

        if (rank == 0) then
                print *, ""
                print *, "Benchmark completed!"
        end if

        call MPI_Finalize(ierr)

        contains

        ! ── CG solver: only mat-vec and dot products are parallelised ──────────
        function cg_mpi(A_loc, b, n, rs, re, rank, nprocs) result(x)
                use mpi
                integer,  intent(in) :: n, rs, re, rank, nprocs
                real,     intent(in) :: A_loc(re-rs+1, n)
                real,     intent(in) :: b(n)

                real    :: alpha, beta, tol
                real    :: rr, rr_new, pAp
                real    :: rr_global, pAp_global
                real    :: r(n), r_new(n), p(n), x(n)
                real    :: Ap(n), Ap_local(re-rs+1)
                integer :: it, ierr

                tol = 1.0e-6
                x   = 0.0
                r   = b
                p   = r
                rr  = dot_product(r, r)

                do it = 1, n

                        ! ── Distributed mat-vec: each rank computes its rows ──
                        Ap_local = matmul(A_loc, p)

                        ! Gather full Ap on all ranks
                        call MPI_Allgather(Ap_local, re-rs+1, MPI_REAL, &
                                           Ap,       re-rs+1, MPI_REAL, &
                                           MPI_COMM_WORLD, ierr)

                        ! ── Dot products with global reduction ────────────────
                        pAp_global = dot_product(p, Ap)   ! same on all ranks (full vectors)

                        alpha = rr / pAp_global

                        x     = x + alpha * p
                        r_new = r - alpha * Ap

                        rr_new = dot_product(r_new, r_new)

                        if (sqrt(rr_new) < tol) then
                                if (rank == 0) print *, "Converged at iteration:", it
                                exit
                        end if

                        beta = rr_new / rr
                        p    = r_new + beta * p
                        r    = r_new
                        rr   = rr_new
                end do
        end function cg_mpi

end program main
