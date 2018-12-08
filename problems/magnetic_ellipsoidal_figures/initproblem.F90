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

  real   ::  a1, a2, a3, alpha, beta, bg_dens

  namelist /PROBLEM_CONTROL/ a1, a2, a3, alpha, beta, bg_dens

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

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

        a1 = rbuff(1)
        a2 = rbuff(2)
        a3 = rbuff(3)
        alpha = rbuff(4)
        beta  = rbuff(5)
        bg_dens = rbuff(6)

      endif
      
  end subroutine read_problem_par
!-----------------------------------------------------------------------------------------------------------------
  subroutine problem_initial_conditions
    
    implicit none

  end subroutine problem_initial_conditions
!----------------------------------------------------------------------------------------------------------------
end module initproblem
