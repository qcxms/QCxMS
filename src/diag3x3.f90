! This is free and unencumbered software released into the public domain.
!
! Anyone is free to copy, modify, publish, use, compile, sell, or
! distribute this software, either in source code form or as a compiled
! binary, for any purpose, commercial or non-commercial, and by any
! means.
!
! In jurisdictions that recognize copyright laws, the author or authors
! of this software dedicate any and all copyright interest in the
! software to the public domain. We make this dedication for the benefit
! of the public at large and to the detriment of our heirs and
! successors. We intend this dedication to be an overt act of
! relinquishment in perpetuity of all present and future rights to this
! software under copyright law.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
! OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
! ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
! OTHER DEALINGS IN THE SOFTWARE.
!
! For more information, please refer to <https://unlicense.org>

!> Analytical implementation of the 3x3 eigenproblem and its solution.
module diag3x3
   implicit none
   private

   public :: eigval3x3, eigvec3x3


   !> Floating point precision
   integer, parameter :: wp = selected_real_kind(15)

   !> 2π/3
   real(wp), parameter :: twothirdpi = 8.0_wp * atan(1.0_wp) / 3.0_wp

   !> Numerical precision
   real(wp), parameter :: eps = epsilon(1.0_wp)


contains


!> Calculates eigenvalues based on the trigonometric solution of A = pB + qI
pure subroutine eigval3x3(a, w)

   !> The symmetric input matrix
   real(wp), intent(in) :: a(3, 3)

   !> Contains eigenvalues on exit
   real(wp), intent(out) :: w(3)

   real(wp) :: q, p, r

   r = a(1, 2) * a(1, 2) + a(1, 3) * a(1, 3) + a(2, 3) * a(2, 3)
   q = (a(1, 1) + a(2, 2) + a(3, 3)) / 3.0_wp
   w(1) = a(1, 1) - q
   w(2) = a(2, 2) - q
   w(3) = a(3, 3) - q
   p = sqrt((w(1) * w(1) + w(2) * w(2) + w(3) * w(3) + 2*r) / 6.0_wp)
   r = (w(1) * (w(2) * w(3) - a(2, 3) * a(2, 3)) &
      & - a(1, 2) * (a(1, 2) * w(3) - a(2, 3) * a(1, 3)) &
      & + a(1, 3) * (a(1, 2) * a(2, 3) - w(2) * a(1, 3))) / (p*p*p) * 0.5_wp

   if (r <= -1.0_wp) then
      r = 0.5_wp * twothirdpi
   else if (r >= 1.0_wp) then
      r = 0.0_wp
   else
      r = acos(r) / 3.0_wp
   end if

   w(3) = q + 2 * p * cos(r)
   w(1) = q + 2 * p * cos(r + twothirdpi)
   w(2) = 3 * q - w(1) - w(3)

end subroutine eigval3x3


!> Calculates eigenvector using an analytical method based on vector cross
!  products.
pure subroutine eigvec3x3(a, w, q)

   !> The symmetric input matrix, destroyed while solving
   real(wp), intent(inout) :: a(3,3)

   !> Contains eigenvalues on exit
   real(wp), intent(out) :: w(3)

   !> Contains eigenvectors on exit
   real(wp), intent(out) :: q(3,3)

   !> Local variables
   real(wp) :: norm, n1, n2, n3, precon
   integer :: i

   w(1) = max(abs(a(1, 1)), abs(a(1, 2)))
   w(2) = max(abs(a(1, 3)), abs(a(2, 2)))
   w(3) = max(abs(a(2, 3)), abs(a(3, 3)))
   precon = max(w(1), max(w(2), w(3)))

   ! null matrix
   if (precon < eps) then
      w(1) = 0.0_wp
      w(2) = 0.0_wp
      w(3) = 0.0_wp
      q(1, 1) = 1.0_wp
      q(2, 2) = 1.0_wp
      q(3, 3) = 1.0_wp
      q(1, 2) = 0.0_wp
      q(1, 3) = 0.0_wp
      q(2, 3) = 0.0_wp
      q(2, 1) = 0.0_wp
      q(3, 1) = 0.0_wp
      q(3, 2) = 0.0_wp
      return
   end if

   norm = 1.0_wp / precon

   a(1, 1) = a(1, 1) * norm
   a(1, 2) = a(1, 2) * norm
   a(2, 2) = a(2, 2) * norm
   a(1, 3) = a(1, 3) * norm
   a(2, 3) = a(2, 3) * norm
   a(3, 3) = a(3, 3) * norm

   ! Calculate eigenvalues
   call eigval3x3(a, w)

   ! Compute first eigenvector
   a(1, 1) = a(1, 1) - w(1)
   a(2, 2) = a(2, 2) - w(1)
   a(3, 3) = a(3, 3) - w(1)

   q(1, 1) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)
   q(2, 1) = a(1, 3) * a(1, 2) - a(1, 1) * a(2, 3)
   q(3, 1) = a(1, 1) * a(2, 2) - a(1, 2) * a(1, 2)
   q(1, 2) = a(1, 2) * a(3, 3) - a(1, 3) * a(2, 3)
   q(2, 2) = a(1, 3) * a(1, 3) - a(1, 1) * a(3, 3)
   q(3, 2) = a(1, 1) * a(2, 3) - a(1, 2) * a(1, 3)
   q(1, 3) = a(2, 2) * a(3, 3) - a(2, 3) * a(2, 3)
   q(2, 3) = a(2, 3) * a(1, 3) - a(1, 2) * a(3, 3)
   q(3, 3) = a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3)
   n1 = q(1, 1) * q(1, 1) + q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)
   n2 = q(1, 2) * q(1, 2) + q(2, 2) * q(2, 2) + q(3, 2) * q(3, 2)
   n3 = q(1, 3) * q(1, 3) + q(2, 3) * q(2, 3) + q(3, 3) * q(3, 3)

   norm = n1
   i = 1
   if (n2 > norm) then
      i = 2
      norm = n1
   end if
   if (n3 > norm) then
      i = 3
   end if

   if (i == 1) then
      norm = sqrt(1.0_wp / n1)
      q(1, 1) = q(1, 1) * norm
      q(2, 1) = q(2, 1) * norm
      q(3, 1) = q(3, 1) * norm
   else if (i == 2) then
      norm = sqrt(1.0_wp / n2)
      q(1, 1) = q(1, 2) * norm
      q(2, 1) = q(2, 2) * norm
      q(3, 1) = q(3, 2) * norm
   else
      norm = sqrt(1.0_wp / n3)
      q(1, 1) = q(1, 3) * norm
      q(2, 1) = q(2, 3) * norm
      q(3, 1) = q(3, 3) * norm
   end if

   ! Robustly compute a right-hand orthonormal set (ev1, u, v)
   if (abs(q(1, 1)) > abs(q(2, 1))) then
      norm = sqrt(1.0_wp / (q(1, 1) * q(1, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = -q(3, 1) * norm
      q(2, 2) = 0.0_wp
      q(3, 2) = +q(1, 1) * norm
   else
      norm = sqrt(1.0_wp / (q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = 0.0_wp
      q(2, 2) = +q(3, 1) * norm
      q(3, 2) = -q(2, 1) * norm
   end if
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   ! Reset A
   a(1, 1) = a(1, 1) + w(1)
   a(2, 2) = a(2, 2) + w(1)
   a(3, 3) = a(3, 3) + w(1)

   ! A*U
   n1 = a(1, 1) * q(1, 2) + a(1, 2) * q(2, 2) + a(1, 3) * q(3, 2)
   n2 = a(1, 2) * q(1, 2) + a(2, 2) * q(2, 2) + a(2, 3) * q(3, 2)
   n3 = a(1, 3) * q(1, 2) + a(2, 3) * q(2, 2) + a(3, 3) * q(3, 2)

   ! A*V, note out of order computation
   a(3, 3) = a(1, 3) * q(1, 3) + a(2, 3) * q(2, 3) + a(3, 3) * q(3, 3)
   a(1, 3) = a(1, 1) * q(1, 3) + a(1, 2) * q(2, 3) + a(1, 3) * q(3, 3)
   a(2, 3) = a(1, 2) * q(1, 3) + a(2, 2) * q(2, 3) + a(2, 3) * q(3, 3)

   ! UT*(A*U) - l2*E
   n1 = q(1, 2) * n1 + q(2, 2) * n2 + q(3, 2) * n3 - w(2)
   ! UT*(A*V)
   n2 = q(1, 2) * a(1, 3) + q(2, 2) * a(2, 3) + q(3, 2) * a(3, 3)
   ! VT*(A*V) - l2*E
   n3 = q(1, 3) * a(1, 3) + q(2, 3) * a(2, 3) + q(3, 3) * a(3, 3) - w(2)

   if (abs(n1) >= abs(n3)) then
      norm = max(abs(n1), abs(n2))
      if (norm > eps) then
         if (abs(n1) >= abs(n2)) then
            n2 = n2 / n1
            n1 = sqrt(1.0_wp / (1.0_wp + n2 * n2))
            n2 = n2 * n1
         else
            n1 = n1 / n2
            n2 = sqrt(1.0_wp / (1.0_wp + n1 * n1))
            n1 = n1 * n2
         end if
         q(1, 2) = n2 * q(1, 2) - n1 * q(1, 3)
         q(2, 2) = n2 * q(2, 2) - n1 * q(2, 3)
         q(3, 2) = n2 * q(3, 2) - n1 * q(3, 3)
      end if
   else
      norm = max(abs(n3), abs(n2))
      if (norm > eps) then
         if (abs(n3) >= abs(n2)) then
            n2 = n2 / n3
            n3 = sqrt(1.0_wp / (1.0_wp + n2 * n2))
            n2 = n2 * n3
         else
            n3 = n3 / n2
            n2 = sqrt(1.0_wp / (1.0_wp + n3 * n3))
            n3 = n3 * n2
         end if
         q(1, 2) = n3 * q(1, 2) - n2 * q(1, 3)
         q(2, 2) = n3 * q(2, 2) - n2 * q(2, 3)
         q(3, 2) = n3 * q(3, 2) - n2 * q(3, 3)
      end if
   end if

   ! Calculate third eigenvector from cross product
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   w(1) = w(1) * precon
   w(2) = w(2) * precon
   w(3) = w(3) * precon

end subroutine eigvec3x3


end module diag3x3

