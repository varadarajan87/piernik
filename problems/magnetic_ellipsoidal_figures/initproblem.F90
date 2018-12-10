!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"

!> \brief Implementation of magnetic ellipsoidal figures of equilibrium.

module initproblem

  implicit none

  private
  public :: read_problem_par, problem_pointers, problem_initial_conditions

  real   ::  a1, a2, a3, alpha, beta, bg_dens, dens_uni

  namelist /PROBLEM_CONTROL/ a1, a2, a3, alpha, beta, bg_dens, dens_uni

contains
!-----------------------------------------------------------------------------------------------------------------
  subroutine problem_pointers

    implicit none

  end subroutine problem_pointers
!-----------------------------------------------------------------------------------------------------------------
  subroutine read_problem_par

    use dataio_pub, only: nh
    use mpisetup,   only: rbuff, master, slave, PIERNIK_MPI_Bcast
    
    implicit none

    a1 = 1.0
    a2 = 1.0
    a3 = 0.6
    alpha = 1.0
    beta  = 1.0
    bg_dens = 1.e-6
    dens_uni = 1.0

    if(master) then
    
        if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = a1
         rbuff(2) = a2
         rbuff(3) = a3
         rbuff(4) = alpha
         rbuff(5) = beta
         rbuff(6) = bg_dens
         rbuff(7) = dens_uni

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

        a1 = rbuff(1)
        a2 = rbuff(2)
        a3 = rbuff(3)
        alpha = rbuff(4)
        beta  = rbuff(5)
        bg_dens = rbuff(6)
        dens_uni = rbuff(7)

      endif
      
  end subroutine read_problem_par
!-----------------------------------------------------------------------------------------------------------------
  subroutine problem_initial_conditions
    
    use cg_leaves,   only: leaves
    use cg_list,     only: cg_list_element
    use constants,   only: xdim, ydim, zdim, one, zero, two
    use grid_cont,   only: grid_container
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use func,        only: ekin, emag
    use units,       only: newtong

    implicit none

    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    class(component_fluid), pointer :: fl

    integer                         :: i,j,k
    real                            :: pres_star, bg_pres
    real                            :: zeta0, Omega, Omega2
    real                            :: C_p, C_B
    real                            :: Z
    real                            :: I, A1, A2, A3, eccty

    eccty = sqrt(one - (a3/a1)**two)
    A1    = (sqrt(one - eccty**two)/e**3.0)*asin(eccty) - (one - eccty**two)/e**two
    A2    = A1
    A3    = two/eccty**two - two*(sqrt(one - eccty**two)/e**3.0)*asin(e)
    I     = a1**two * A1 + a2**two * A2 + a3**two * A3 
    C_p   = (-(I - A1*a1**two - A2*a2**two - A3*a3**two)) + &
               zeta0**2*(alpha**two * a1**two + beta**two * a2**two)/(two*(alpha + beta)**two) - &
                ( (zeta0**two/(two*(alpha + beta))) + one )*(alpha*a1**two + beta*a2**two)
    !C_p = ( A2 - zeta0**two* ( ((a1**two * a2**two)/(two*(a1**two+a2**two)**two))  + Z*(a1/a2)**two  ) )*a2**two - ( I )
    C_B  = 4.0*pi*dens_uni*(alpha*a1**two + beta*a2**two)
    !Determine after fixing units
    Omega2 = two*sqrt((one - eccty**two))*(3.0 - two*eccty**two)*asin(eccty)/e**3.0 - 6.0*(one - e**two)/e**two ! In units of pi.G.\rho
    Omega  = sqrt(Omega2) ! In units of sqrt(pi.G.\rho)
    zeta0  = two*Omega
    Z      = alpha/zeta0**two

    ! Background is basically a vacuum. For numerical considerations we prescribe a 
    ! small density with no initial motion and no magentic fields.
    ! We do not consider the background fluid to interact with the ellipsoid at t = 0. 
    ! Hence from hydrostatic equilibrium the pressure is also uniform. 
    bg_pres = bg_dens

    fl  => flind%ion
    cgl => leaves%first

    do while (associated(cgl))
       
       cg => cgl%cg

       do k = cg%ks, cg%ke
          do j = cg%js, cg%je
             do i = cg%is, cg%ie
                
                cg%u(fl%imx,i,j,k) = zero
                cg%u(fl%imy,i,j,k) = zero
                cg%u(fl%imz,i,j,k) = zero
                cg%b(xdim,i,j,k)   = zero
                cg%b(ydim,i,j,k)   = zero
                cg%b(zdim,i,j,k)   = zero

                ! Pressure of the star. Eq.(15) 
                pres_star = dens_uni*( (I - A1*cg%x(i)*cg*x(i) - A2*cg%y(j)*cg%y(j) -A3*cg%z(k)*cg%z(k)) - &
                                            zeta0**two*(alpha**two * cg%x(i)*cg%x(i) + beta**two * cg%y(j)*cg%y(j))/(two*(alpha + beta)**two) + &
                                            ( (zeta0**2/(2*(alpha + beta))) + one )*(alpha*cg%x(i)*cg%x(i) + beta*cg%y(j)*cg%y(j)) + C_p )
                
                if(pres_star .gt. bg_pres) then
                   
                   cg%u(fl%idn,i,j,k) = dens_uni
                   cg%u(fl%imx,i,j,k) = -cg%u(fl%idn,i,j,k)*zeta0*(beta*cg%y(j))/(alpha + beta)
                   cg%u(fl%imy,i,j,k) = cg%u(fl%idn,i,j,k)*zeta0*(alpha*cg%x(i))/(alpha + beta)
                   cg%b(zdim,i,j,k)   = sqrt(two*(C_B - 4.0*pi*dens_uni*(alpha*cg%x(i)*cg%x(i) + beta*cg%y(j)*cg%y(j))))
                   cg%u(fl%ien,i,j,k) = pres_star/fl%gam_1 + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                        emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                   
                else 
                   
                   cg%u(fl%idn,i,j,k) = bg_dens
                   cg%u(fl%ien,i,j,k) = bg_pres/fl%gam_1 + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                        emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                   
                   
                end if
                
             end do
          end do
       end do
       cgl => cgl%nxt
    end do

  end subroutine problem_initial_conditions
!----------------------------------------------------------------------------------------------------------------
end module initproblem
