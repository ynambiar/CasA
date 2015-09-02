
;+
;rd_casa.pro
;	read in all the files created from CasA data reduction
;
;vinay k (2013-sep-20)
;	streamlined (2013-nov-18)
;-

;	initialize
rootdir='/data/darmok/CasA/'
;
obsids=[198,1046,1952,5196,6744,6745,9117,10935,10936,14229]
obsids=[198,1952,5196,6744,6745,9117,10935,10936,14229]	;get rid of 1046 because it's a grating obs
obsids=[198,1952,5196,6745,9117,10935,10936,14229]	;get rid of 6744 because it's ACIS-I
nobs=n_elements(obsids) & sobsids=string(obsids,'(i5.5)') & cobsids=strtrim(obsids,2)
;
bands=['soft','medium','hard','broad'] & nbands=n_elements(bands)
;
binsiz=[4,8,16,32,64] & nbinsz=n_elements(binsiz) & sbinsz=string(binsiz,'(i3.3)')
numx=1024/binsiz & numy=768/binsiz

;	skip if already read in
if keyword_set(already_read_in) then goto,skiprd

;	define the outputs to hold the arrays read in
for i=0L,nbinsz-1L do begin
  bb=sbinsz[i]
  jnk=execute("ctimg"+bb+"=lonarr(nobs,nbands,numx[i],numy[i])")
  jnk=execute("emapimg"+bb+"=fltarr(nobs,nbands,numx[i],numy[i])")
  jnk=execute("flximg"+bb+"=fltarr(nobs,nbands,numx[i],numy[i])")
endfor
obsdate=strarr(nobs) & exptim=dblarr(nobs)

;	read in the files for each ObsID
for i=0L,nobs-1L do begin		;{for each ObsID
  for j=0L,nbinsz-1L do begin		;{for each bin size
    bb=sbinsz[j]
    for k=0L,nbands-1L do begin		;{for each band

      ctfil=rootdir+'/'+cobsids[i]+'/repro/f'+sbinsz[j]+'/'+bands[k]+'.img'
      emapfil=rootdir+'/'+cobsids[i]+'/repro/f'+sbinsz[j]+'/'+bands[k]+'.expmap'
      flxfil=rootdir+'/'+cobsids[i]+'/repro/f'+sbinsz[j]+'/'+bands[k]+'_flux.img'
      tmp=file_search(ctfil,count=nctfil) & if nctfil eq 0 then message,ctfil+': not found'
      tmp=file_search(emapfil,count=nemapfil) & if nemapfil eq 0 then message,emapfil+': not found'
      tmp=file_search(flxfil,count=nflxfil) & if nflxfil eq 0 then message,flxfil+': not found'

      spawn,'ls -l '+ctfil+' '+emapfil+' '+flxfil

      ctimg=mrdfits(ctfil,0,hct)
      emapimg=mrdfits(emapfil,0,hemap)
      flximg=mrdfits(flxfil,0,hflx)

      if j eq 0 and k eq 0 then begin
        obsdate[i]=sxpar(hct,'DATE-OBS') & exptim[i]=sxpar(hct,'EXPOSURE')
        print,'ObsID: ',cobsids[i]
        print,'Observation Date: ',obsdate[i]
        print,'Exposure [s]: ',exptim[i]
      endif

      jnk=execute("ctimg"+bb+"[i,k,*,*]=ctimg[0:numx[j]-1,0:numy[j]-1]")
      jnk=execute("emapimg"+bb+"[i,k,*,*]=emapimg[0:numx[j]-1,0:numy[j]-1]")
      jnk=execute("flximg"+bb+"[i,k,*,*]=flximg[0:numx[j]-1,0:numy[j]-1]")

      if !d.name eq 'X' then begin
        window,0,xsize=numx[0],ysize=numy[0] & tvscl,alog10(ctimg>0.9)
        window,1,xsize=numx[0],ysize=numy[0] & tvscl,emapimg
        window,2,xsize=numx[0],ysize=numy[0] & tvscl,alog10(flximg>(1./exptim[i]/300.))
      endif

    endfor				;K=0,NBANDS-1}
  endfor				;J=0,NBINSZ-1}
endfor					;I=0,NOBS-1}
save,file='rd_casa.save'
skiprd: already_read_in=1
spawn,'ls -l rd_casa.save'

end
