;+
;pp_casa.pro
;	post-process the files read in from CasA reduction
;
;vinay k (2013-sep-22)
;	debugged and enhanced (2013-nov18)
;-

;	initialize
if n_elements(flximg004) eq 0 then begin
  message,'reading in output of Cas A reduction',/informational
  restore,'/data/darmok/CasA/rd_casa.save'
endif
;
peasecolr & loadct,9 & peasecolr

;	compute some simple metrics
for ibin=0L,nbinsz-1L do begin
  bb=sbinsz[ibin]	;this is the binning scale

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("ctimg=ctimg"+bb+" & emapimg=emapimg"+bb+" & flximg=flximg"+bb)

  ;	first, compute some useful stuff on pixel-wise basis
  ; sum across ObsIDs to get total counts in each pixel
  sumimg=lonarr(nbands,numx[ibin],numy[ibin]) & for i=0,nbands-1 do sumimg[i,*,*]=total(reform(ctimg[*,i,*,*]),1)
  ; average flux across ObsIDs to get flux in each pixel
  avgimg=fltarr(nbands,numx[ibin],numy[ibin]) & for i=0,nbands-1 do avgimg[i,*,*]=total(reform(flximg[*,i,*,*]),1)/nobs

  ; make a mask based on which pixels cover parts of the remnant
  ; a pixel must have at least 1 count in each pixel in all bands in all observations (otherwise colors will have infinities),
  ; must have 200 counts over all observations, and must have an average flux >5e-6 ph/s/cm^2 at bin=4
  ;	to get to 5e-6, do: plot,histogram(alog10(avgimg[3,*,*]),min=-7,max=-4,bin=0.1)
  ;	and there is a minimum at bin~15, so push to 17 to minimize background contamination, which corresponds to 5e-6
  minct=intarr(numx[ibin],numy[ibin])+max(ctimg) & for i1=0L,nobs-1L do for i2=0,nbands-2 do minct=minct < (ctimg[i1,i2,*,*])
  oy=where(minct gt 0 and reform(sumimg[3,*,*] gt 200 and reform(avgimg[3,*,*]) gt 5e-6*(binsiz[ibin]/binsiz[0])^2),moy,complement=xoy,ncomplement=mxoy)
  maskimg=bytarr(numx[ibin],numy[ibin])+1 & if mxoy gt 0 then maskimg[xoy]=0
  	;	OY (inside) and XOY (outside) are special arrays -- DO NOT OVERWRITE BELOW!
  tvscl,maskimg	;<- this just plots the mask

  ;	this part is not used in determining the result, but is left here as potential future extension
  ; compute stddev on fluxes at each pixel
  sdfimg=fltarr(nbands,numx[ibin],numy[ibin]) & for i=0,nbands-1 do for ix=0L,numx[ibin]-1 do for jy=0L,numy[ibin]-1 do $
  	sdfimg[i,*,*]=stddev(flximg[*,i,ix,jy])
  ; compute statistical error on flux at each point
  errimg=sqrt(ctimg) & o0=where(ctimg gt 0,mo0) & if mo0 gt 0 then errimg[o0]=(errimg[o0]/float(ctimg[o0]))*flximg[o0]
  ; compute a chisq value for flux variations for each pixel
  chiimg=fltarr(nbands,numx[ibin],numy[ibin])
  for i=0,nbands-1 do begin
    tmp=fltarr(numx[ibin],numy[ibin])
    for ix=0L,numx[ibin]-1 do begin
      for jy=0L,numy[ibin]-1 do begin
        tmp1=reform(flximg[*,i,ix,jy]) & tmp2=reform(errimg[*,i,ix,jy]) & tmp3=avgimg[i,ix,jy]
        tmp[ix,jy]=total((tmp1-tmp3)^2/tmp2^2)
      endfor
    endfor
    if mxoy gt 0 then tmp[xoy]=0.
    chiimg[i,*,*]=tmp
  endfor

  ;	now compute some useful stuff on observation-wise basis, suitably masked for interesting pixels
  ; total counts in each ObsID
  sumobs=lonarr(nobs,nbands) & for i=0L,nobs-1L do for j=0L,nbands-1L do sumobs[i,j]=total(reform(ctimg[i,j,*,*]*maskimg))
  ; total flux in each ObsID
  flxobs=fltarr(nobs,nbands) & for i=0L,nobs-1L do for j=0L,nbands-1L do flxobs[i,j]=total(reform(flximg[i,j,*,*]*maskimg))
  ; stddev of flux in each ObsID
  sdfobs=fltarr(nobs,nbands) & for i=0L,nobs-1L do for j=0L,nbands-1L do begin
    tmp=reform(flximg[i,j,*,*]) & if moy gt 0 then sdfobs[i,j]=stddev(tmp[oy])
  endfor

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("sumimg"+bb+"=sumimg & avgimg"+bb+"=avgimg & oy"+bb+"=oy & xoy"+bb+"=xoy & minct"+bb+"=minct & maskimg"+bb+"=maskimg & sdfimg"+bb+"=sdfimg & errimg"+bb+"=errimg & chiimg"+bb+"=chiimg & sumobs"+bb+"=sumobs & flxobs"+bb+"=flxobs & sdfobs"+bb+"=sdfobs")

  ;	color images
  ;	make flux images for each band, suitably masked
  flxS=reform(flximg[*,0,*,*]) & flxM=reform(flximg[*,1,*,*]) & flxH=reform(flximg[*,2,*,*])
  for j=0L,nobs-1L do flxS[j,*,*]=flxS[j,*,*]*maskimg
  for j=0L,nobs-1L do flxM[j,*,*]=flxM[j,*,*]*maskimg
  for j=0L,nobs-1L do flxH[j,*,*]=flxH[j,*,*]*maskimg
  ;csm=alog10(float(flxS)/float(flxM)) ;& oo=where(finite(csm) eq 0,moo) & if moo gt 0 then csm[oo]=-30.	;these will have NaNs
  ;cmh=alog10(float(flxM)/float(flxH)) ;& oo=where(finite(cmh) eq 0,moo) & if moo gt 0 then cmh[oo]=-30.	;these will have NaNs
  csm=fltarr(nobs,numx[ibin],numy[ibin]) & cmh=csm
  for i=0L,nobs-1L do begin
    tmp=fltarr(numx[ibin],numy[ibin]) & tmp1=reform(flxS[i,*,*]) & tmp2=reform(flxM[i,*,*]) & if moy gt 0 then tmp[oy]=alog10(tmp1[oy]/tmp2[oy]) & csm[i,*,*]=tmp
    tmp=fltarr(numx[ibin],numy[ibin]) & tmp1=reform(flxM[i,*,*]) & tmp2=reform(flxH[i,*,*]) & if moy gt 0 then tmp[oy]=alog10(tmp1[oy]/tmp2[oy]) & cmh[i,*,*]=tmp
    	;NaN free, but that means must use OY and XOY
  endfor

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("csm"+bb+"=csm & cmh"+bb+"=cmh")

  ;	now analyze the ratios to find abnormal excursions in ratio values, for both ratios

  ; initialize for ratio 1, Soft/Medium
  z1log=csm

  ; average per pixel
  z1logavpix=fltarr(numx[ibin],numy[ibin]) &  for ix=0L,numx[ibin]-1L do for jy=0L,numy[ibin]-1L do $
  	z1logavpix[ix,jy]=mean(z1log[*,ix,jy],/nan)

  ; stddev per pixel
  z1logsdpix=fltarr(numx[ibin],numy[ibin]) & for ix=0L,numx[ibin]-1L do for jy=0L,numy[ibin]-1L do $
  	if total(finite(z1log[*,ix,jy]) ne 0) gt 1 then z1logsdpix[ix,jy]=stddev(z1log[*,ix,jy],/nan)

  ;; chisq per pixel (not used)
  ;z1logchipix=fltarr(numx[ibin],numy[ibin]) & for ix=0L,numx[ibin]-1L do for jy=0L,numy[ibin]-1L do $
  ;	if total(finite(z1log[*,ix,jy]) ne 0) gt 1 then z1logchipix[ix,iy]=total((z1log[*,ix,jy]-z1logavpix[ix,jy])^2)

  ; average per observation and stddev of values per observation
  z1logavobs=fltarr(nobs) & z1logsdobs=fltarr(nobs)
  for j=0L,nobs-1L do begin
    tmp=reform(z1log[j,*,*])
    if moy gt 0 then z1logavobs[j]=mean(tmp[oy],/nan)
    if moy gt 0 then z1logsdobs[j]=stddev(tmp[oy],/nan)
  endfor

  ; renormalize to get light curve for each pixel by removing average
  z1norm=z1log & for j=0L,nobs-1L do z1norm[j,*,*]=z1norm[j,*,*]-z1logavpix[*,*]	;this renormalizes all ratios to "unity" so large fluctuations can be detected as N stddev deviations

  ; stddev for each observation for renormalized values
  z1normsdobs=fltarr(nobs)
  for j=0L,nobs-1L do begin
    tmp=reform(z1norm[j,*,*]) & if moy gt 0 then z1normsdobs[j]=stddev(tmp[oy],/nan)
  endfor

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("z1logavpix"+bb+"=z1logavpix & z1logsdpix"+bb+"=z1logsdpix & z1logavobs"+bb+"=z1logavobs & z1logsdobs"+bb+"=z1logsdobs & z1norm"+bb+"=z1norm & z1normsdobs"+bb+"=z1normsdobs")

  ; renormalize to recenter the ratio for each observation
  z1renorm=z1log & for j=0L,nobs-1L do z1renorm[j,*,*]=z1renorm[j,*,*]-z1logavobs[j]	;this moves the ratios to the "middle" for each ObsID
  z1sdevobs=fltarr(nobs,numx[ibin],numy[ibin])
  ; define the thresholds to call something variable
  z1thrp=3*z1logsdobs & z1thrm=-3*z1logsdobs
  ; bitmap image, for each observation, of all pixels which show values very different from the rest of the pixels
  z1bitsobs=intarr(nobs,numx[ibin],numy[ibin])
  for j=0L,nobs-1L do begin
    tmp=reform(z1renorm[j,*,*]) & oo=where(tmp gt z1thrp[j] or tmp lt z1thrm[j],moo,complement=xoo,ncomplement=mxoo)
    if moo gt 0 then tmp[oo]=1
    if mxoo gt 0 then tmp[xoo]=0
    if mxoy gt 0 then tmp[xoy]=0
    z1bitsobs[j,*,*]=tmp
    z1sdevobs[j,*,*]=(z1renorm[j,*,*]/z1thrp[j]) > (z1renorm[j,*,*]/z1thrm[j])
  endfor

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("z1renorm"+bb+"=z1renorm & z1thrp"+bb+"=z1thrp & z1thrm"+bb+"=z1thrm & z1bitsobs"+bb+"=z1bitsobs & z1sdevobs"+bb+"=z1sdevobs")

  ; for each lightcurve at each pixel, we now have an estimate of its scatter in z1logsdpix
  ; now pick out the extreme values of the scatter, and that shows which pixels are the most variable
  z1bitsp=intarr(numx[ibin],numy[ibin])
  z1sdevpix=fltarr(numx[ibin],numy[ibin])
  tmp=z1logsdpix[oy] & mtmp=mean(tmp) & stmp=stddev(tmp)
  oho=where(z1logsdpix gt mtmp+3*stmp,moho)
  if moho gt 0 then z1bitsp[oho]=1
  if moho gt 0 then z1sdevpix[oho]=z1logsdpix[oho]/stmp

  ; summed bitmaps, for pixel-wise and observation-wise interesting pixels
  z1bitso=total(z1bitsobs,1)

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("z1bitsp"+bb+"=z1bitsp & z1bitso"+bb+"=z1bitso & z1sdevpix"+bb+"=z1sdevpix")
  ;	z1bitsp contains pixels flagged for large variability across all observations
  ;	z1bitso contains pixels flagged for abnormal values within one observation
  ;	z1sdevpix contains ratio of lc stddev relative to threshold stddev

  ;------------------------------------------------------------------------
  ; for Medium/Hard, same as Soft/Medium
  z2log=cmh

  z2logavpix=fltarr(numx[ibin],numy[ibin]) &  for ix=0L,numx[ibin]-1L do for jy=0L,numy[ibin]-1L do $
  	z2logavpix[ix,jy]=mean(z2log[*,ix,jy],/nan)

  z2logsdpix=fltarr(numx[ibin],numy[ibin]) & for ix=0L,numx[ibin]-1L do for jy=0L,numy[ibin]-1L do $
  	if total(finite(z2log[*,ix,jy]) ne 0) gt 1 then z2logsdpix[ix,jy]=stddev(z2log[*,ix,jy],/nan)

  z2logavobs=fltarr(nobs) & z2logsdobs=fltarr(nobs)

  for j=0L,nobs-1L do begin
    tmp=reform(z2log[j,*,*])
    if moy gt 0 then z2logavobs[j]=mean(tmp[oy],/nan)
    if moy gt 0 then z2logsdobs[j]=stddev(tmp[oy],/nan)
  endfor

  z2norm=z2log & for j=0L,nobs-1L do z2norm[j,*,*]=z2norm[j,*,*]-z2logavpix[*,*]	;this renormalizes all ratios to "unity" so large fluctuations can be detected as N stddev deviations
  z2normsdobs=fltarr(nobs)
  for j=0L,nobs-1L do begin
    tmp=reform(z2norm[j,*,*]) & if moy gt 0 then z2normsdobs[j]=stddev(tmp[oy],/nan)
  endfor

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("z2logavpix"+bb+"=z2logavpix & z2logsdpix"+bb+"=z2logsdpix & z2logavobs"+bb+"=z2logavobs & z2logsdobs"+bb+"=z2logsdobs & z2norm"+bb+"=z2norm & z2normsdobs"+bb+"=z2normsdobs")

  z2renorm=z2log & for j=0L,nobs-1L do z2renorm[j,*,*]=z2renorm[j,*,*]-z2logavobs[j]	;this moves the ratios to the "middle" for each ObsID
  z2thrp=3*z2logsdobs & z2thrm=-3*z2logsdobs
  z2sdevobs=fltarr(nobs,numx[ibin],numy[ibin])
  z2bitsobs=intarr(nobs,numx[ibin],numy[ibin])
  for j=0L,nobs-1L do begin
    tmp=reform(z2renorm[j,*,*]) & oo=where(tmp gt z2thrp[j] or tmp lt z2thrm[j],moo,complement=xoo,ncomplement=mxoo)
    if moo gt 0 then tmp[oo]=1
    if mxoo gt 0 then tmp[xoo]=0
    if mxoy gt 0 then tmp[xoy]=0
    z2bitsobs[j,*,*]=tmp
    z2sdevobs[j,*,*]=(z2renorm[j,*,*]/z2thrp[j]) > (z2renorm[j,*,*]/z2thrm[j])
  endfor

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("z2renorm"+bb+"=z2renorm & z2thrp"+bb+"=z2thrp & z2thrm"+bb+"=z2thrm & z2bitsobs"+bb+"=z2bitsobs & z2sdevobs"+bb+"=z2sdevobs")

  ; for each lightcurve at each pixel, we now have an estimate of its scatter in z1logsdpix
  ; now pick out the extreme values of the scatter, and that shows which pixels are the most variable
  z2bitsp=intarr(numx[ibin],numy[ibin])
  z2sdevpix=fltarr(numx[ibin],numy[ibin])
  tmp=z2logsdpix[oy] & mtmp=mean(tmp) & stmp=stddev(tmp)
  oho=where(z2logsdpix gt mtmp+3*stmp,moho)
  if moho gt 0 then z2bitsp[oho]=1
  if moho gt 0 then z2sdevpix[oho]=z2logsdpix[oho]/stmp

  ; summed bitmaps, for pixel-wise and observation-wise interesting pixels
  z2bitso=total(z2bitsobs,1)

  ;	the next line just stores the arrays with names tagged by the scale size so that they may be recalled at the end
  jnk=execute("z2bitsp"+bb+"=z2bitsp & z2bitso"+bb+"=z2bitso & z2sdevpix"+bb+"=z2sdevpix")
  ;	z2bitsp contains pixels flagged for large variability across all observations
  ;	z2bitso contains pixels flagged for abnormal values within one observation

  if !d.name eq 'X' then begin
    window,1 & tvscl,rebin(z1bitso gt 0,numx[0]*2,numy[0]*2) & print,total(z1bitso)
    window,2 & tvscl,rebin(z2bitsp gt 0,numx[0]*2,numy[0]*2) & print,total(z1bitsp)
  endif

endfor

;--------------------------------------------------------------------------------
;	write out all the identified pixels
;	note: these are in image coordinates.
;	can convert to WCS by loading them in over corresponding bin-size image in ds9,
;	then selecting all the regions and saving them in ciao/wcs format.

openw,u004o,'casa004_idxo_1.reg',/get_lun
openw,u004os,'casa004_idxo_sig_1.txt',/get_lun
o1=where(Z1BITSO004 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSO004,o1)
for i=0,mo1-1 do printf,u004o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u004os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z1SDEVOBS004[*,ii[0,i],ii[1,i]]),2)
close,u004o & free_lun,u004o & close,u004os & free_lun,u004os
openw,u004p,'casa004_idxp_1.reg',/get_lun
openw,u004ps,'casa004_idxp_sig_1.txt',/get_lun
o1=where(Z1BITSP004 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSP004,o1)
for i=0,mo1-1 do printf,u004p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u004ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z1SDEVPIX004[ii[0,i],ii[1,i]],2)
close,u004p & free_lun,u004p & close,u004ps & free_lun,u004ps

openw,u008o,'casa008_idxo_1.reg',/get_lun
openw,u008os,'casa008_idxo_sig_1.txt',/get_lun
o1=where(Z1BITSO008 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSO008,o1)
for i=0,mo1-1 do printf,u008o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u008os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z1SDEVOBS008[*,ii[0,i],ii[1,i]]),2)
close,u008o & free_lun,u008o & close,u008os & free_lun,u008os
openw,u008p,'casa008_idxp_1.reg',/get_lun
openw,u008ps,'casa008_idxp_sig_1.txt',/get_lun
o1=where(Z1BITSP008 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSP008,o1)
for i=0,mo1-1 do printf,u008p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u008ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z1SDEVPIX008[ii[0,i],ii[1,i]],2)
close,u008p & free_lun,u008p & close,u008ps & free_lun,u008ps

openw,u016o,'casa016_idxo_1.reg',/get_lun
openw,u016os,'casa016_idxo_sig_1.txt',/get_lun
o1=where(Z1BITSO016 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSO016,o1)
for i=0,mo1-1 do printf,u016o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u016os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z1SDEVOBS016[*,ii[0,i],ii[1,i]]),2)
close,u016o & free_lun,u016o & close,u016os & free_lun,u016os
openw,u016p,'casa016_idxp_1.reg',/get_lun
openw,u016ps,'casa016_idxp_sig_1.txt',/get_lun
o1=where(Z1BITSP016 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSP016,o1)
for i=0,mo1-1 do printf,u016p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u016ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z1SDEVPIX016[ii[0,i],ii[1,i]],2)
close,u016p & free_lun,u016p & close,u016ps & free_lun,u016ps

openw,u032o,'casa032_idxo_1.reg',/get_lun
openw,u032os,'casa032_idxo_sig_1.txt',/get_lun
o1=where(Z1BITSO032 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSO032,o1)
for i=0,mo1-1 do printf,u032o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u032os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z1SDEVOBS032[*,ii[0,i],ii[1,i]]),2)
close,u032o & free_lun,u032o & close,u032os & free_lun,u032os
openw,u032p,'casa032_idxp_1.reg',/get_lun
openw,u032ps,'casa032_idxp_sig_1.txt',/get_lun
o1=where(Z1BITSP032 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSP032,o1)
for i=0,mo1-1 do printf,u032p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u032ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z1SDEVPIX032[ii[0,i],ii[1,i]],2)
close,u032p & free_lun,u032p & close,u032ps & free_lun,u032ps

openw,u064o,'casa064_idxo_1.reg',/get_lun
openw,u064os,'casa064_idxo_sig_1.txt',/get_lun
o1=where(Z1BITSO064 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSO064,o1)
for i=0,mo1-1 do printf,u064o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u064os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z1SDEVOBS064[*,ii[0,i],ii[1,i]]),2)
close,u064o & free_lun,u064o & close,u064os & free_lun,u064os
openw,u064p,'casa064_idxp_1.reg',/get_lun
openw,u064ps,'casa064_idxp_sig_1.txt',/get_lun
o1=where(Z1BITSP064 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z1BITSP064,o1)
for i=0,mo1-1 do printf,u064p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u064ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z1SDEVPIX064[ii[0,i],ii[1,i]],2)
close,u064p & free_lun,u064p & close,u064ps & free_lun,u064ps

;------------------------------------------------------------------------

openw,u004o,'casa004_idxo_2.reg',/get_lun
openw,u004os,'casa004_idxo_sig_2.txt',/get_lun
o1=where(Z2BITSO004 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSO004,o1)
for i=0,mo1-1 do printf,u004o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u004os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z2SDEVOBS004[*,ii[0,i],ii[1,i]]),2)
close,u004o & free_lun,u004o & close,u004os & free_lun,u004os
openw,u004p,'casa004_idxp_2.reg',/get_lun
openw,u004ps,'casa004_idxp_sig_2.txt',/get_lun
o1=where(Z2BITSP004 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSP004,o1)
for i=0,mo1-1 do printf,u004p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u004ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z2SDEVPIX004[ii[0,i],ii[1,i]],2)
close,u004p & free_lun,u004p & close,u004ps & free_lun,u004ps

openw,u008o,'casa008_idxo_2.reg',/get_lun
openw,u008os,'casa008_idxo_sig_2.txt',/get_lun
o1=where(Z2BITSO008 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSO008,o1)
for i=0,mo1-1 do printf,u008o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u008os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z2SDEVOBS008[*,ii[0,i],ii[1,i]]),2)
close,u008o & free_lun,u008o & close,u008os & free_lun,u008os
openw,u008p,'casa008_idxp_2.reg',/get_lun
openw,u008ps,'casa008_idxp_sig_2.txt',/get_lun
o1=where(Z2BITSP008 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSP008,o1)
for i=0,mo1-1 do printf,u008p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u008ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z2SDEVPIX008[ii[0,i],ii[1,i]],2)
close,u008p & free_lun,u008p & close,u008ps & free_lun,u008ps

openw,u016o,'casa016_idxo_2.reg',/get_lun
openw,u016os,'casa016_idxo_sig_2.txt',/get_lun
o1=where(Z2BITSO016 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSO016,o1)
for i=0,mo1-1 do printf,u016o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u016os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z2SDEVOBS016[*,ii[0,i],ii[1,i]]),2)
close,u016o & free_lun,u016o & close,u016os & free_lun,u016os
openw,u016p,'casa016_idxp_2.reg',/get_lun
openw,u016ps,'casa016_idxp_sig_2.txt',/get_lun
o1=where(Z2BITSP016 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSP016,o1)
for i=0,mo1-1 do printf,u016p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u016ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z2SDEVPIX016[ii[0,i],ii[1,i]],2)
close,u016p & free_lun,u016p & close,u016ps & free_lun,u016ps

openw,u032o,'casa032_idxo_2.reg',/get_lun
openw,u032os,'casa032_idxo_sig_2.txt',/get_lun
o1=where(Z2BITSO032 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSO032,o1)
for i=0,mo1-1 do printf,u032o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u032os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z2SDEVOBS032[*,ii[0,i],ii[1,i]]),2)
close,u032o & free_lun,u032o & close,u032os & free_lun,u032os
openw,u032p,'casa032_idxp_2.reg',/get_lun
openw,u032ps,'casa032_idxp_sig_2.txt',/get_lun
o1=where(Z2BITSP032 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSP032,o1)
for i=0,mo1-1 do printf,u032p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u032ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z2SDEVPIX032[ii[0,i],ii[1,i]],2)
close,u032p & free_lun,u032p & close,u032ps & free_lun,u032ps

openw,u064o,'casa064_idxo_2.reg',/get_lun
openw,u064os,'casa064_idxo_sig_2.txt',/get_lun
o1=where(Z2BITSO064 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSO064,o1)
for i=0,mo1-1 do printf,u064o,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u064os,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(max(Z2SDEVOBS064[*,ii[0,i],ii[1,i]]),2)
close,u064o & free_lun,u064o & close,u064os & free_lun,u064os
openw,u064p,'casa064_idxp_2.reg',/get_lun
openw,u064ps,'casa064_idxp_sig_2.txt',/get_lun
o1=where(Z2BITSP064 gt 0,mo1) & if mo1 gt 0 then ii=array_indices(Z2BITSP064,o1)
for i=0,mo1-1 do printf,u064p,'box('+strtrim(ii[0,i],2)+','+strtrim(ii[1,i],2)+',1.0,1.0,0)'
for i=0,mo1-1 do printf,u064ps,strtrim(ii[0,i],2)+'	'+strtrim(ii[1,i],2)+'	'+strtrim(Z2SDEVPIX064[ii[0,i],ii[1,i]],2)
close,u064p & free_lun,u064p & close,u064ps & free_lun,u064ps

save,file='pp_casa.save'
spawn,'ls -l pp_casa.save casa???_idx?_?.reg casa???_idx?_sig_?.txt'

end
