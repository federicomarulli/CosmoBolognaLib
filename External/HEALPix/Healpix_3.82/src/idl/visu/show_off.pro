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
pro show_off

flip=0
ps= 1

map = findgen(48)
triangle=create_struct('coord','G','ra',[0,80,0],'dec',[40,45,65])

cartview,map,/online,res=45,gra=[45,30],rot=[10,20,30],$
pysize=250,outline=triangle, $
title='Cartesian cylindrical (full sky)',subtitle='cartview'  ,$
    ,ps=ps,/pre

gnomview,map,/online,res=25,gra=[45,30],rot=[10,20,30],$
outline=triangle, $
title='Gnomic projection',subtitle='gnomview'     ,$
    ,ps=ps,/pre

mollview,map,/online,gra=[45,30],rot=[10,20,30],$
outline=triangle, $
title='Mollweide projection',subtitle='mollview'      ,$
      ,ps=ps,/pre

orthview,map,/online,gra=[45,30],rot=[10,20,30],$
outline=triangle, $
title='Orthographic projection',subtitle='orthview'    ,$
     ,ps=ps,/pre

; 

; orthview,'test.fits',pol=3,rot=[80,15,30],/grat,units=' ',
; title='Temperature+Polarization',/pre,ps='/tmp/orthpol.ps'

return
end
