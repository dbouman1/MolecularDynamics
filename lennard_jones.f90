! lennard_jones function calculates the forces between particles and total potential energy

! lennard_jones(r_min, eps, L, r_c, num_bins, p, N, d, force, Ep, vir, c_bin)

! input:
!   r_min       -> lj parameter
!   eps         -> lj parameter
!   L           -> size of the box
!   r_c         -> cutoff length for interaction
!   num_bins    -> number of bins for correlation function
!   p           -> position components matrix
!   N           -> number of particles (automatic)
!   d           -> dimensionality (automatic)

! output:
!   force   -> force components matrix
!   Ep      -> potential energy
!   vir     -> factor for the virial theorem pressure calculation
!   c_bin   -> bin for the correlation

subroutine lennard_jones(r_min, eps, L, r_c, num_bins, p, N, d, force, Ep, vir, c_bin)
    implicit none

    integer, intent(in) :: num_bins 
    integer, intent(inout) :: N, d;
    real(8), intent(in) :: r_min, eps, r_c, L
    real(8), intent(in), dimension(N, d) :: p

    real(8), intent(out), dimension(N, d) :: force
    integer, intent (out), dimension(num_bins) :: c_bin
    real(8), intent(out) :: Ep, vir

    real(8), dimension(d) :: r
    real(8) :: r_abs, force_fact
    integer :: i, j, bin_num

    do i=1,N
        do j=1,N
            r = p(i,:) - p(j,:) - L*ANINT((p(i,:) - p(j,:))/L)
			r_abs = sqrt(sum(r**2))

            if ((r_abs .NE. 0) .and. (r_abs < r_c))  then

                ! bin the distances
                bin_num = CEILING(num_bins*r_abs/r_c)
                c_bin(bin_num) = c_bin(bin_num) + 1

                ! potential energy
                Ep = Ep + &
                    0.5*eps*((r_min**12/r_abs**12) - 2.0*(r_min**6/r_abs**6))
                
                ! force matrix
				force_fact = 12.0*eps*((r_min**12/r_abs**14) - (r_min**6/r_abs**8))
                force(i,:) = force(i,:) + force_fact*r
					
				! virial factor
				vir = vir - 0.5*force_fact*r_abs**2
            end if
        end do
    end do

end subroutine lennard_jones
