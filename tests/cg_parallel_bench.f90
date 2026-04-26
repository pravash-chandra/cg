program main 
        !$use omp_lib
        integer, dimension(5) :: n_values
        integer :: n, test_idx
        real, allocatable, dimension(:,:) :: A 
        real, allocatable, dimension(:) :: b, x
        integer :: i, j
        real :: start_time, end_time, elapsed_time
        
        ! Create an array of size 5 ranging from 5000 to 25000 for testing the speed
        n_values = [5000, 10000, 15000, 20000, 25000]
        
        ! Print header
        print *, "====== CG Solver Benchmarking (OpenMP Parallel) ======"
        print *, ""
        print '(A15, A15)', "Matrix Size", "Time (sec)"
        print *, "-------------------------------------------"
        
        ! Loop through each test size
        do test_idx = 1, 5
                n = n_values(test_idx)
                
                ! Allocate arrays
                if (allocated(A)) deallocate(A)
                if (allocated(b)) deallocate(b)
                if (allocated(x)) deallocate(x)
                allocate(A(n,n))
                allocate(b(n))
                allocate(x(n))
        
                ! Dense symmetric positive definite matrix for benchmarking
                do i = 1, n
                    do j = 1, n
                        A(i,j) = 1.0 / (abs(i-j) + 1.0)
                    end do
                end do
                do i = 1, n
                    A(i,i) = A(i,i) + real(n)
                end do
        
                ! Create b vector of same size
                b = 1.0
                
                ! Call CG solver with timing
                call cpu_time(start_time)
                x = cg(A,b,n)
                call cpu_time(end_time)
                
                elapsed_time = end_time - start_time
                print '(I15, F15.6)', n, elapsed_time
        
        end do
        
        print *, ""
        print *, "Benchmark completed!"

        contains
        function cg(A, b,n) result(x)
                ! inputs declarations 
                integer , intent(in) :: n
                real , dimension(n,n), intent(in) :: A
                real , dimension(n) , intent(in) :: b

                real :: alpha , beta
                real, dimension(n) :: r , r_new, p, Ap
                real, dimension(n) :: x
                integer :: i, j
                real :: tolerance, rr, pAp

                ! initialize x to 0 vector of size n
                tolerance = 1.0e-6
                x = 0.0
                r = b
                p = r
                
                ! since in CG in krylov space it is guranted to find the solution 
                ! before no of iterations = size of  matrix
                do i = 1,n 
                        ! Matrix-vector product: Ap = A * p (parallelized)
                        !$omp parallel do private(j)
                        do j = 1, n
                                Ap(j) = dot_prod(A(j,:), p)
                        end do
                        !$omp end parallel do
                        
                        ! Dot products
                        rr = dot_prod(r, r)
                        pAp = dot_prod(p, Ap)
                        
                        alpha = rr / pAp
                        
                        ! Update x: x = x + alpha * p (parallelized)
                        !$omp parallel do
                        do j = 1, n
                                x(j) = x(j) + alpha * p(j)
                        end do
                        !$omp end parallel do
                        
                        ! Update r: r_new = r - alpha * Ap (parallelized)
                        !$omp parallel do
                        do j = 1, n
                                r_new(j) = r(j) - alpha * Ap(j)
                        end do
                        !$omp end parallel do
                        
                        if (sqrt(dot_prod(r_new, r_new)) < tolerance) then
                                print *, "Converged at iteration:", i
                                exit
                        end if
                        
                        beta = dot_prod(r_new, r_new) / rr
                        
                        ! Update p: p = r_new + beta * p (parallelized)
                        !$omp parallel do
                        do j = 1, n
                                p(j) = r_new(j) + beta * p(j)
                        end do
                        !$omp end parallel do
                        
                        r = r_new
                end do
        end function cg 
        
        function dot_prod(a, b) result(res)
                real, dimension(:), intent(in) :: a, b
                real :: res
                res = sum(a * b)
        end function dot_prod

end program main
