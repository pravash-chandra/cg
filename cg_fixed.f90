program main 
        integer :: n
        real, allocatable, dimension(:,:) :: A 
        real, allocatable, dimension(:) :: b, x
        integer :: i, j
        
        ! Create A a diagonal matrix so that it is guaranteed to be symmetric 
        ! positive definite

        n = 5
        ! Allocate arrays
        allocate(A(n,n))
        allocate(b(n))
        allocate(x(n))
        
        ! diagonal matrix for testing if the solution is correct and match against 
        ! analytical solutions
        A = 0.0
        do i = 1, n
                A(i,i) = real(i) * 2.0  ! Diagonal entries: 2, 4, 6, 8, 10
        end do

        ! ! Dense symmetric positive definite matrix for benchmarking
        ! do i = 1, n
        !     do j = 1, n
        !         A(i,j) = 1.0 / (abs(i-j) + 1.0)
        !     end do
        ! end do
        ! A(i,i) = A(i,i) + real(n)

        
        ! Create b vector of same size
        b = 1.0
        
        ! Call CG solver
        x = cg(A,b,n)
        
        ! Print the result in formatted output
        print *, "Solution vector x:"
        print '(5F12.6)', x
        print *, ""
        print *, "Matrix A:"
        do i = 1, n
                print '(5F12.6)', A(i,:)
        end do

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
                        alpha = dot_product(r, r) / dot_product(p, matmul(A, p))
                        x = x + alpha * p
                        r_new = r - alpha * matmul(A, p)
                        
                        if (sqrt(dot_product(r_new, r_new)) < tolerance) then
                                print *, "Converged at iteration:", i
                                exit
                        end if
                        
                        beta = dot_product(r_new, r_new) / dot_product(r, r)
                        p = r_new + beta * p
                        r = r_new
                end do
        end function cg 
end program main
