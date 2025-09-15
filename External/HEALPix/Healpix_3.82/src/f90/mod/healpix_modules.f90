!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2013 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.sourceforge.net
!
!-----------------------------------------------------------------------------
module healpix_modules
! Healpix meta module, giving access to all specific modules


  use healpix_types
  use extension
  use misc_utils
  use utilities
  use paramfile_io
  use long_intrinsic

  use healpix_fft
  use num_rec
!  use m_indmed
  use statistics
  use rngmod
  use coord_v_convert

  use bit_manipulation
  use pix_tools
  use head_fits
  use fitstools
  use alm_tools
!  use mpi_alm_tools
  use udgrade_nr
  use mask_tools

  use ran_tools
  use obsolete

  public

end module healpix_modules
