program main 
        integer, dimension(5) :: n_values
        integer :: n, test_idx
        real, allocatable, dimension(:,:) :: A 
        real, allocatable, dimension(:) :: b, x
        integer :: i, j
        real :: start_time, end_time, elapsed_time
        
        ! Create an array of size 5 ranging from 5000 to 25000 for testing the speed
        n_values = [5000, 10000, 15000, 20000, 25000]
        
        ! Print header
        print *, "====== CG Solver Benchmarking (Sequential) ======"
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
                real, dimension(n) :: r , r_new, p
                real, dimension(n) :: x
                integer :: i
                real :: tolerance

                ! initialize x to 0 vector of size n
                tolerance = 1.0e-6
                x = 0.0
                r = b
                p = r
                
                ! since in CG in krylov space it is guranted to find the solution 
                ! before no of iterations = size of  matrix
                do i = 1,n 
                        alpha = dot_prod(r, r) / dot_prod(p, mat_prod(A, p))
                        x = x + alpha * p
                        r_new = r - alpha * mat_prod(A, p)
                        
                        if (sqrt(dot_prod(r_new, r_new)) < tolerance) then
                                print *, "Converged at iteration:", i
                                exit
                        end if
                        
                        beta = dot_prod(r_new, r_new) / dot_prod(r, r)
                        p = r_new + beta * p
                        r = r_new
                end do
        end function cg 
        
        function dot_prod(a, b) result(res)
                real, dimension(:), intent(in) :: a, b
                real :: res
                res = sum(a * b)
        end function dot_prod
        
        function mat_prod(A, v) result(result_vec)
                real, dimension(:,:), intent(in) :: A
                real, dimension(:), intent(in) :: v
                real, dimension(size(v)) :: result_vec
                integer :: i, j, n
                
                n = size(v)
                do i = 1, n
                        result_vec(i) = 0.0
                        do j = 1, n
                                result_vec(i) = result_vec(i) + A(i,j) * v(j)
                        end do
                end do
        end function mat_prod

end program main
