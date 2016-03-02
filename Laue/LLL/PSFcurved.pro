pro psf

totpix=float(1)
totpix = 1024.0 ;squared detector with 1024 x 1024 pixels
centerx =  abs(totpix/2) 
centery =  abs(totpix/2)
peak_refl = 0.14 ; peak reflectivity value
incoming_source = 400.   ;  photons/(cm^2 s)
focal = 20. ; (m)
hc=12.39842
dhkl=3.135531
EB=150.0



ringradius = (hc * focal) / (dhkl * EB ) * 100; (m)
pixel_dimension = 0.3 ; (mm)
external_curvature = focal * 2; (m)
internal_curvature = 2.56 * external_curvature ; (m) 
dim1 = 11.7  ; radial (mm)
dim2 = 15.0  ; tangential (mm)
dim3 = 5.0   ; thickness  (mm)
ncristals = 12; 10*ringradius*2*3.1417/dim2 ; total number of crystals in that ring
crys_surface = dim1 * dim2 / 100 ; cm^2 total surface of the crystal 
ph = incoming_source * crys_surface

print, ph
tB=asin(hc / (2. * dhkl * EB))

alpha=2*atan((dim1/2)/(external_curvature * 1000.))    
Emin = hc / (2 * dhkl * sin(tB + alpha/2))  
Emax = hc / (2 * dhkl * sin(tB - alpha/2))  
deltaE = Emax - Emin


; calcolo le dimensioni sul detector
det1 = dim2 
det2 = dim3/4

; trasformo in pixel
px1 = floor(det1/pixel_dimension)
px2 = floor(det2/pixel_dimension)



startx = centerx - px1/2
stopx = centerx + px1/2

starty = centery - px2/2
stopy = centery + px2/2

cris_image=fltarr(totpix,totpix)


for i=startx,stopx do for j=starty,stopy do cris_image[i,j] = ph * peak_refl * (deltaE) / px1*px2; (det1*det2) 

; noise creation
detector_noise = 0.4 * RANDOMN(POISSON,[totpix,totpix])
for i = 0, totpix-1 DO FOR j = 0, totpix-1 DO BEGIN 
if float(j-centerx)^2 + float(i-centery)^2 gt ((totpix/2)^2) then begin
   detector_noise[j,i]=-5000
endif 

endfor

signal =  cris_image
for n = 1, ncristals-1 DO BEGIN 
cris=fltarr(totpix,totpix)
cris = rot(cris_image,n*360/ncristals)
signal =  signal + cris
endfor

randnumb = 2 * (RANDOMU(S, totpix, totpix) - 0.5)

ph_noise = randnumb * sqrt(abs(signal))

signalplusnoise = signal + ph_noise

a=signalplusnoise ; congrid(signalplusnoise,totpix,totpix)  

a_zoom=rebin(a(450:577,450:577), 512,512,/sample)

;#############################################################
;thisDevice = !D.NAME
;SET_PLOT, 'Z'
;LOADCT, 5
;TVLCT, red, green, blue, /GET
;SET_PLOT, thisDevice
;thisImage = BYTSCL(a_zoom)
;s = SIZE(thisImage)
;image3d = BYTARR(3, s(1), s(2))
;image3d(0, *, *) = red(thisImage)
;image3d(1, *, *) = green(thisImage)
;image3d(2, *, *) = blue(thisImage)
;WRITE_JPEG, 'myimage.jpeg', image3d, TRUE=1, QUALITY=100
;##############################################################
;TVLCT, r,g,b  
;TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
;TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
;TVLCT, INDGEN(256), INDGEN(256), INDGEN(256)
;TVLCT, [[255], [255], [0]], 100
LOADCT, 13
SET_PLOT, 'ps', /interpolate
LOADCT, 13
DEVICE, decomposed=0, FILE='24bit.ps', BITS=8 , /color
shade_surf, a_zoom, /device, ax=45.,az=30., shades=bytscl(a_zoom)



;dataColors = !D.TABLE_SIZE-2
;redColor = !D.TABLE_SIZE-2
;greenColor = !D.TABLE_SIZE-1
;LOADCT, 13;, NCOLORS=datacolors
;TVLCT, [255, 0], [0, 255], [0,0], !D.TABLE_SIZE-4 ; Red and Green.
;TVLCT, red, green, blue, /GET
;POLYFILL, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=0
;!P.COLOR = 0
;!P.BACKGROUND = 255 
;DEVICE, FILE='24bit.ps', BITS=8
;shade_surf, a_zoom, /device, ax=45.,az=30., shades=bytscl(a_zoom), /color
;isurface, a_zoom;, /device, ax=45.,az=30., shades=bytscl(a_zoom), /color, /skirt
;TV, [[[2]], [[3]], [[1]]], TRUE=3
;thisDevice = !D.NAME
;SET_PLOT, thisDevice
;device,/close
;set_plot,'x'

;!P.POSITION=[0.1,0.15,0.8,0.85]
;shade_surf,a_zoom,ax=90,az=0,charsize=2., title='Curved crystal PSF'; ,  COLOR = 300
;shade_surf, a_zoom, ax=45.,az=45., shade=bytscl(alog(a_zoom+0.001))
;surface,a_zoom,xstyle=4,ystyle=4,zstyle=4,/noerase
;show3,a

;#####################################################################
;window,0, retain=2, xsize=512,ysize=512
;device, decomposed=0
;loadct, 13
;shade_surf, a_zoom, ax=45.,az=45., shade=bytscl(alog(a_zoom+0.001))
;#####################################################################

;TV, a_zoom, True=1
;entry_device = !d.name
;set_plot, 'PS'
;device, /color,  filename='FullScatterPlot.ps'
;!P.Multi=0
;shade_surf, a_zoom, ax=45.,az=45., shade=bytscl(alog(a_zoom+0.001))
;device,/close
;set_plot,'x'


;set_plot,'ps'
;device, filename='postscriptfilename.ps', /COLOR, BITS=8
;shade_surf, a_zoom, ax=45.,az=45., shade=bytscl(alog(a_zoom+0.001))
;device,/close
;set_plot,'x'


;TVLCT, red, green, blue, /GET
;TVLCT, 5, 5, 5, 0
;tvscl, (a_zoom)



;cgSurface, aa


end
