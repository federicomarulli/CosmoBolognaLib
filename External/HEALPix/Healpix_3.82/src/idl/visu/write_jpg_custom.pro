pro write_jpg_custom, file, image, red, green, blue, $
                  order = order, $
                  progressive = progressive, $
                  quality = quality, $
                  true = true, $
                  unit = unit
;+
; write_jpg_custom, file, image [, red, green, blue, ORDER=, PROGRESSIVE=,
; QUALITY=, TRUE=, UNIT=]
;
;  wrapper around write_jpeg with syntax closer to write_png
;
; 2012-01-20
;-


sz =size(image)
if (sz[0] eq 3) then begin
    write_jpeg, file, image, $
                  order = order, $
                  progressive = progressive, $
                  quality = quality, $
                  true = true, $
                  unit = unit
endif else begin
    image3d = bytarr(3,sz[1],sz[2], /nozero)
    image3d[0, *, *] = red  [image]
    image3d[1, *, *] = green[image]
    image3d[2, *, *] = blue [image]
    write_jpeg, file, image3d, $
                  order = order, $
                  progressive = progressive, $
                  quality = quality, $
                  true = 1, $
                  unit = unit
endelse

return
end
