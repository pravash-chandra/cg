program cg_openmp_benchmark
        integer :: n, i, j, test_idx
        real, allocatable, dimension(:,:) :: A 
        real, allocatable, dimension(:) :: b, x
        real :: start_time, end_time, elapsed_time
        integer :: sizes(6)
        integer :: num_threads
        
        ! Matrix sizes to test: 1000 to 15000
        sizes = [1000, 2500, 5000, 7500, 10000, 15000]
        
        ! Get number of threads
        !$omp parallel
        !$omp single
        num_threads = omp_get_num_threads()
        !$omp end single
        !$omp end parallel
        
        print *, "=========================================================================="
        print '(A, I3, A)', "OPENMP CG BENCHMARK - Dense Matrix (", num_threads, " threads)"
        print *, "=========================================================================="
        print *, ""
        print '(A15, A20, A15, A15)', "Matrix Size", "Time (seconds)", "Iterations", "Status"
        print *, "=========================================================================="
        
        do test_idx = 1, 6
                n = sizes(test_idx)
                
                ! Allocate arrays
                allocate(A(n,n))
                allocate(b(n))
                allocate(x(n))
                
                ! Initialize b
                b = 1.0
                
                ! Create dense symmetric positive definite matrix
                ! Using parallelized initialization
                A = 0.0
                !$omp parallel do private(j)
                do i = 1, n
                        do j = 1, n
                                A(i,j) = 1.0 / (abs(i-j) + 1.0)
                        end do
                end do
                !$omp end parallel do
                
                ! Make diagonal dominant for stability
                !$omp parallel do
                do i = 1, n
                        A(i,i) = A(i,i) + real(n) * 0.1
                end do
                !$omp end parallel do
                
                ! Start timing
                call cpu_time(start_time)
                
                ! Call CG solver
                call cg_solver(A, b, x, n)
                
                ! End timing
                call cpu_time(end_time)
                elapsed_time = end_time - start_time
                
                ! Print results
                print '(I15, F20.6, A15, A15)', n, elapsed_time, "", "✓ Complete"
                
                ! Cleanup
                deallocate(A, b, x)
        end do
        
        print *, "=========================================================================="
        print *, ""
        print '(A, I3, A)', "OpenMP Configuration: ", num_threads, " threads"
        print *, "Notes:"
        print *, "  - Dense symmetric positive definite matrix"
        print *, "  - Convergence tolerance: 1.0e-6"
        print *, "  - Times measured using CPU time"
        print *, "  - Set OMP_NUM_THREADS to change thread count"
        print *, "=========================================================================="

        contains
        subroutine cg_solver(A, b, x, n)
                integer, intent(in) :: n
                real, dimension(n,n), intent(in) :: A
                real, dimension(n), intent(in) :: b
                real, dimension(n), intent(out) :: x
                
                real :: alpha, beta
                real, dimension(n) :: r, r_new, p, Ap
                integer :: i, k
                real :: tolerance, rr, pAp
                
                ! Initialize
                tolerance = 1.0e-6
                x = 0.0
                r = b
                p = r
                
                ! CG iterations
                do i = 1, n
                        ! Matrix-vector product (parallelized)
                        !$omp parallel do private(k)
                        do k = 1, n
                                Ap(k) = dot_product(A(k,:), p)
                        end do
                        !$omp end parallel do
                        
                        rr = dot_product(r, r)
                        pAp = dot_product(p, Ap)
                        
                        if (abs(pAp) < 1.0e-14) exit
                        
                        alpha = rr / pAp
                        
                        ! Update x (parallelized)
                        !$omp parallel do
                        do k = 1, n
                                x(k) = x(k) + alpha * p(k)
                        end do
                        !$omp end parallel do
                        
                        ! Update r (parallelized)
                        !$omp parallel do
                        do k = 1, n
                                r_new(k) = r(k) - alpha * Ap(k)
                        end do
                        !$omp end parallel do
                        
                        if (sqrt(dot_product(r_new, r_new)) < tolerance) then
                                exit
                        end if
                        
                        beta = dot_product(r_new, r_new) / rr
                        
                        ! Update p (parallelized)
                        !$omp parallel do
                        do k = 1, n
                                p(k) = r_new(k) + beta * p(k)
                        end do
                        !$omp end parallel do
                        
                        r = r_new
                end do
        end subroutine cg_solver
        
end program cg_openmp_benchmark
