; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2013  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.sourceforge.net
;
; -----------------------------------------------------------------------------
function sub_compute_vertices, z, z_nv, z_sv, phi, phi_nv, phi_sv, hdelta_phi
;+
; function called by pix2vec_ring (pix2vec_nest) to compute the
; 3D vectors pointing toward the vertices of a pixel from their angles 
;
;
;-

np = n_elements(z)
vertex = dblarr(np, 3, 4)

sth = SQRT((1.0d0-z)*(1.0d0+z))

sth_nv = sqrt((1.0d0-z_nv)*(1.0d0+z_nv))
vertex[*,*,0] = [[sth_nv*cos(phi_nv)],[sth_nv*sin(phi_nv)],[z_nv]] ; north vertex

phi_wv = phi - hdelta_phi
vertex[*,*,1] = [[sth*cos(phi_wv)],[sth*sin(phi_wv)],[z]] ; west vertex

sth_sv = sqrt((1.0d0-z_sv)*(1.0d0+z_sv))
vertex[*,*,2] = [[sth_sv*cos(phi_sv)],[sth_sv*sin(phi_sv)],[z_sv]] ; south vertex

phi_ev = phi + hdelta_phi
vertex[*,*,3] = [[sth*cos(phi_ev)],[sth*sin(phi_ev)],[z]] ; east vertex





return, vertex
end
