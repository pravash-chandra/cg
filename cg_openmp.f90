program main 
        integer                        :: n
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
        
        ! Initialize A to diagonal matrix with positive entries
        A = 0.0
        do i = 1, n
                A(i,i) = real(i) * 2.0  ! Diagonal entries
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
        function cg(A, b, n) result(x)
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
                
                ! CG iterations
                do i = 1, n 
                        ! Matrix-vector product: Ap = A * p (parallelized)
                        !$omp parallel do private(j)
                        do j = 1, n
                                Ap(j) = dot_product(A(j,:), p)
                        end do
                        !$omp end parallel do
                        
                        ! Dot products (not parallelized - small overhead)
                        rr = dot_product(r, r)
                        pAp = dot_product(p, Ap)
                        
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
                        
                        ! Check convergence
                        if (sqrt(dot_product(r_new, r_new)) < tolerance) then
                                print *, "Converged at iteration:", i
                                exit
                        end if
                        
                        ! Beta calculation
                        beta = dot_product(r_new, r_new) / rr
                        
                        ! Update p: p = r_new + beta * p (parallelized)
                        !$omp parallel do
                        do j = 1, n
                                p(j) = r_new(j) + beta * p(j)
                        end do
                        !$omp end parallel do
                        
                        r = r_new
                end do
        end function cg 
end program main
