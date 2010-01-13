!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module lsvar 


contains
  ! 
  ! ----------------------------------------------------------------------
  ! subroutine varlscale : computes the time variation of the sst and
  ! the divergence for astex case

  subroutine varlscale(time_in,case_name,sst,div,u0,v0)

  use grid, only : umean, vmean,nzp 
  use forc, only :  t_ls,div_ls,sst_ls,ugeo_ls,vgeo_ls

  implicit none

   real, optional, intent (in) :: time_in
   real, optional, intent (inout) :: sst,div
   real, optional,  intent (inout) :: u0(nzp),v0(nzp)
   character (len=5), intent (in) :: case_name

   real                  :: tlag,tfrac,ugeo,vgeo
   real            :: del_sst,del_div,del_ugeo,del_vgeo
   integer         :: tcnt
   integer               :: k


   tlag=(t_ls(2)-t_ls(1))*3600.
!   print *,'tlag', tlag

   tfrac=(time_in/tlag)-int(time_in/tlag)
   tcnt=int(time_in/tlag)+2


   del_sst=(sst_ls(tcnt)-sst_ls(tcnt-1))
   del_div=(div_ls(tcnt)-div_ls(tcnt-1))
   del_ugeo=(ugeo_ls(tcnt)-ugeo_ls(tcnt-1))
   del_vgeo=(vgeo_ls(tcnt)-vgeo_ls(tcnt-1))
   
   sst=sst_ls(tcnt-1)+del_sst*tfrac
   div=div_ls(tcnt-1)+del_div*tfrac
   div=div*1.e-6
   ugeo=ugeo_ls(tcnt-1)+del_ugeo*tfrac
   vgeo=vgeo_ls(tcnt-1)+del_vgeo*tfrac
  
 if (trim(case_name) == 'astex') then
    do k=1,nzp
    u0(k) = ugeo - umean
    v0(k) = vgeo - vmean
    end do
 end if

!if (hour_frac .ge. 20. .and. hour_frac .ge. 20.2 ) then
!print *, 'u0 vo', u0,v0
!end if

   !print *,'lsvar', time_in, tcnt, tfrac, sst, div

   end subroutine varlscale

end module lsvar   


