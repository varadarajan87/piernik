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
!   The literature references are MNRAS 416, L75-L79 (2011), hereafter Kawa11.
!   Ellipsoidal figures of equilibrium, S. Chandrasekhar, hereafter SC.
!<

module initproblem

  implicit none

  private
  public :: read_problem_par, problem_pointers, problem_initial_conditions

  real   ::  a1, a2, a3, alpha, beta, bg_dens, bg_pres, dens_uni

  namelist /PROBLEM_CONTROL/ a1, a2, a3, alpha, beta, bg_dens, bg_pres, dens_uni

contains
!-----------------------------------------------------------------------------------------------------------------
  subroutine problem_pointers

    use gravity,               only: grav_pot_3d

    implicit none

    grav_pot_3d => ellipsoid_grav_pot_3d

  end subroutine problem_pointers
!-----------------------------------------------------------------------------------------------------------------
  subroutine read_problem_par

    use dataio_pub, only: nh
    use mpisetup,   only: rbuff, master, slave, PIERNIK_MPI_Bcast
    
    implicit none

    a1 = 1.e6
    a2 = 1.e6
    a3 = 0.6e6
    alpha = 1.0
    beta  = 1.0
    bg_dens = 1.e-8
    bg_pres = 1.e-8
    dens_uni = 1.e14

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
         rbuff(7) = bg_pres
         rbuff(8) = dens_uni

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

        a1 = rbuff(1)
        a2 = rbuff(2)
        a3 = rbuff(3)
        alpha = rbuff(4)
        beta  = rbuff(5)
        bg_dens = rbuff(6)
        bg_pres = rbuff(7)
        dens_uni = rbuff(8)

      endif
      
  end subroutine read_problem_par
!-----------------------------------------------------------------------------------------------------------------
  subroutine problem_initial_conditions
    
    use cg_leaves,   only: leaves
    use cg_list,     only: cg_list_element
    use constants,   only: xdim, ydim, zdim, one, zero, two, pi
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
    real                            :: pres_star
    real                            :: zeta0, Omega, Omega2
    real                            :: C_p, C_B
    real                            :: Z
    real                            :: I_ellip, AA1, AA2, AA3, eccty

    ! SC Eq(37), Ch 3
    eccty = sqrt(one - (a3/a1)**two) ! eccty = 0.7141428429

    ! SC Eq(36), Ch 3
    AA1    = (sqrt(one - eccty**two)/eccty**3.0)*asin(eccty) - (one - eccty**two)/eccty**two
    AA2    = AA1
    AA3    = two/eccty**two - two*(sqrt(one - eccty**two)/eccty**3.0)*asin(eccty)
    I_ellip     = a1**two * AA1 + a2**two * AA2 + a3**two * AA3 ! Defined after Eq(16) in Kawa11

    ! SC Eq(6), Ch 5
    Omega2 = (two*sqrt((one - eccty**two))*(3.0 - two*eccty**two)*asin(eccty)/eccty**3.0 - 6.0*(one - eccty**two)/eccty**two)*pi*newtong*dens_uni 
    Omega  = sqrt(Omega2) ! T = 2*pi/Omega = 2.55265390e-3 sec, Omega = 2461.432513 sec^{-1}

    ! Kawa11 Eq(18)
    zeta0  = two*Omega
    Z      = alpha/zeta0**two

    ! Kawa11 from Eq(15)
    C_p   = (-pi*newtong*dens_uni*(I_ellip - AA1*a1**two - AA2*a2**two - AA3*a3**two)) + &
               zeta0**2*(alpha**two * a1**two + beta**two * a2**two)/(two*(alpha + beta)**two) - &
                ( (zeta0**two/(two*(alpha + beta))) + one )*(alpha*a1**two + beta*a2**two)
    !C_p = ( AA2 - zeta0**two* ( ((a1**two * a2**two)/(two*(a1**two+a2**two)**two))  + Z*(a1/a2)**two  ) )*a2**two - ( I ) ! Kawa11 Eq(19)
    ! Kawa11 from Eq(12)
    C_B  = 4.0*pi*dens_uni*(alpha*a1**two + beta*a2**two)
   
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
                pres_star = dens_uni*( pi*newtong*dens_uni*(I_ellip - AA1*cg%x(i)*cg%x(i) - AA2*cg%y(j)*cg%y(j) -AA3*cg%z(k)*cg%z(k)) - &
                                            zeta0**two*(alpha**two * cg%x(i)*cg%x(i) + beta**two * cg%y(j)*cg%y(j))/(two*(alpha + beta)**two) + &
                                            ( (zeta0**2/(2*(alpha + beta))) + one )*(alpha*cg%x(i)*cg%x(i) + beta*cg%y(j)*cg%y(j)) + C_p )
                
                if(pres_star .gt. bg_pres) then
                   
                   cg%u(fl%idn,i,j,k) = dens_uni
                   cg%u(fl%imx,i,j,k) = -cg%u(fl%idn,i,j,k)*zeta0*(beta*cg%y(j))/(alpha + beta) ! Kawa11 Eq(6)
                   cg%u(fl%imy,i,j,k) = cg%u(fl%idn,i,j,k)*zeta0*(alpha*cg%x(i))/(alpha + beta) ! Kawa11 Eq(6)
                   cg%b(zdim,i,j,k)   = sqrt(two*(C_B - 4.0*pi*dens_uni*(alpha*cg%x(i)*cg%x(i) + beta*cg%y(j)*cg%y(j)))) ! Kawa11 Eq(12)
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

  subroutine ellipsoid_grav_pot_3d

     use cg_leaves,   only: leaves
     use cg_list,     only: cg_list_element
     use constants,   only: pi, one, two
     use grid_cont,   only: grid_container
     use fluidindex,  only: flind
     use fluidtypes,  only: component_fluid
     use units,       only: newtong

    implicit none

    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    class(component_fluid), pointer :: fl

    integer                         :: i, j, k
    real                            :: I_ellip, AA1, AA2, AA3, eccty

    eccty = sqrt(one - (a3/a1)**two)
    AA1    = (sqrt(one - eccty**two)/eccty**3.0)*asin(eccty) - (one - eccty**two)/eccty**two
    AA2    = AA1
    AA3    = two/eccty**two - two*(sqrt(one - eccty**two)/eccty**3.0)*asin(eccty)
    I_ellip     = a1**two * AA1 + a2**two * AA2 + a3**two * AA3 

    fl  => flind%ion
    cgl => leaves%first

    do while (associated(cgl))
       
       cg => cgl%cg

       if(.not. cg%is_old) then

          do k = cg%ks, cg%ke
             do j = cg%js, cg%je
                do i = cg%is, cg%ie
                   cg%u(fl%idn,i,:,k) = dens_uni
                   cg%gp(i,j,k) = -pi*newtong*cg%u(fl%idn,i,j,k)*(I - AA1*cg%x(i)*cg%x(i) - AA2*cg%y(j)*cg%y(j) - AA3*cg%z(k)*cg%z(k)) ! SC Eq(40) Ch 3
                end do
             end do
          end do

       end if

         cgl => cgl%nxt

    end do
    
  end subroutine ellipsoid_grav_pot_3d

!----------------------------------------------------------------------------------------------------------------

end module initproblem
