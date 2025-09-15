function process_reserved_colors, rgb, color=color

if size(color, /tname) eq 'STRING' then begin
    rgb1 = htmlcolor2rgb(color[0])
endif else begin
    if n_elements(color) eq 3 then begin
        rgb1 = color
    endif else begin
        rgb1 = reform(rgb[color,0:2],3)
    endelse
endelse

return, rgb1
end



pro reserved_colors, red, green, blue, bad_color=bad_color, bg_color=bg_color, fg_color=fg_color, idx_bwg=idx_bwg

routine = 'reserved_colors'
syntax1 = routine+' RGB,     bad_color=, bg_color=, fg_color=, idx_bwg='
syntax2 = routine+' R, G, B, bad_color=, bg_color=, fg_color=, idx_bwg='

np = n_params()
if (np ne 1 && np ne 3) then begin
    print,'expected either 1 or 3 arguments:'
    print,syntax1
    print,syntax2
    return
endif
if np eq 1 then begin
    rgb = red
    red   = rgb[*,0]
    green = rgb[*,1]
    blue  = rgb[*,2]
endif

; set up some specific definitions
; reserve first colors for Black, White and Neutral grey
idx_black = 0B & idx_white = 1B   & idx_grey = 2B   & idx_bwg = [idx_black, idx_white, idx_grey]
col_black = 0B & col_white = 255B & col_grey = 175B & col_bwg = [col_black, col_white, col_grey]
red  [idx_bwg] = col_bwg
green[idx_bwg] = col_bwg
blue [idx_bwg] = col_bwg

rgb = [[red],[green],[blue]]
myrgb = rgb
if defined(bad_color) then myrgb[idx_grey , 0:2]= process_reserved_colors( rgb, color=bad_color)
if defined(bg_color)  then myrgb[idx_white, 0:2]= process_reserved_colors( rgb, color=bg_color )
if defined(fg_color)  then myrgb[idx_black, 0:2]= process_reserved_colors( rgb, color=fg_color )
if np eq 1 then begin
    rgb = myrbg
endif else begin
    red = myrgb[*,0] & green = myrgb[*,1] & blue = myrgb[*,2]
endelse


return
end

