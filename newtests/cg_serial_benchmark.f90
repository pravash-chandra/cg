program cg_serial_benchmark
        integer :: n, i, j, test_idx
        real, allocatable, dimension(:,:) :: A 
        real, allocatable, dimension(:) :: b, x
        real :: start_time, end_time, elapsed_time
        integer :: sizes(6)
        
        ! Matrix sizes to test: 1000 to 15000
        sizes = [1000, 2500, 5000, 7500, 10000, 15000]
        
        print *, "=========================================================================="
        print *, "SERIAL CG BENCHMARK - Dense Matrix (1000 to 15000)"
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
                ! A(i,j) = 1.0 / (abs(i-j) + 1.0)
                ! This creates a realistic sparse-like structure
                A = 0.0
                do i = 1, n
                        do j = 1, n
                                A(i,j) = 1.0 / (abs(i-j) + 1.0)
                        end do
                end do
                
                ! Make diagonal dominant for stability
                do i = 1, n
                        A(i,i) = A(i,i) + real(n) * 0.1
                end do
                
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
        print *, "Notes:"
        print *, "  - Dense symmetric positive definite matrix"
        print *, "  - Convergence tolerance: 1.0e-6"
        print *, "  - Times measured using CPU time"
        print *, "=========================================================================="

        contains
        subroutine cg_solver(A, b, x, n)
                integer, intent(in) :: n
                real, dimension(n,n), intent(in) :: A
                real, dimension(n), intent(in) :: b
                real, dimension(n), intent(out) :: x
                
                real :: alpha, beta
                real, dimension(n) :: r, r_new, p
                integer :: i
                real :: tolerance, rr, pAp
                
                ! Initialize
                tolerance = 1.0e-6
                x = 0.0
                r = b
                p = r
                
                ! CG iterations
                do i = 1, n
                        rr = dot_product(r, r)
                        pAp = dot_product(p, matmul(A, p))
                        
                        if (abs(pAp) < 1.0e-14) exit
                        
                        alpha = rr / pAp
                        x = x + alpha * p
                        r_new = r - alpha * matmul(A, p)
                        
                        if (sqrt(dot_product(r_new, r_new)) < tolerance) then
                                exit
                        end if
                        
                        beta = dot_product(r_new, r_new) / rr
                        p = r_new + beta * p
                        r = r_new
                end do
        end subroutine cg_solver
        
end program cg_serial_benchmark
