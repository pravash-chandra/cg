module cg_solver
    implicit none
    
contains
    
    subroutine solve_cg(A, b, x, n) bind(C, name="solve_cg")
        use iso_c_binding
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: A(n,n)
        real(c_double), intent(in) :: b(n)
        real(c_double), intent(out) :: x(n)
        
        real(c_double) :: alpha, beta
        real(c_double) :: r(n), r_new(n), p(n)
        integer :: i
        real(c_double) :: tolerance
        
        ! Initialize
        tolerance = 1.0d-6
        x = 0.0d0
        r = b
        p = r
        
        ! CG iterations
        do i = 1, n
            alpha = dot_product(r, r) / dot_product(p, matmul(A, p))
            x = x + alpha * p
            r_new = r - alpha * matmul(A, p)
            
            if (sqrt(dot_product(r_new, r_new)) < tolerance) then
                exit
            end if
            
            beta = dot_product(r_new, r_new) / dot_product(r, r)
            p = r_new + beta * p
            r = r_new
        end do
        
    end subroutine solve_cg
    
end module cg_solver
