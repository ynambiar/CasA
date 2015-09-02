
#!/bin/tcsh -f

foreach dir ( `find /data/darmok/CasA -type d -name repro` )

  cd $dir
  set evtfil = `find $dir -name \*_repro_evt2.fits`
  set asollis = `find $dir -name \*asol1.lis`
  set bpixfil = `find $dir -name \*_repro_bpix1.fits`
  set maskfil = `find $dir -name \*msk1.fits`
  set pbkfil = `find $dir -name \*pbk0.fits`

  ls -l $evtfil $asollis $bpixfil $maskfil $pbkfil

  fluximage infile=${evtfil}"[sky=region(/data/darmok/CasA/casa_box.reg)]" outroot=f004/ bands="csc,broad" xygrid="" binsize=4 asolfile=@${asollis} badpixfile=$bpixfil maskfile=$maskfil pbkfile=$pbkfil cleanup=no clobber=yes verbose=2
  fluximage infile=${evtfil}"[sky=region(/data/darmok/CasA/casa_box.reg)]" outroot=f008/ bands="csc,broad" xygrid="" binsize=8 asolfile=@${asollis} badpixfile=$bpixfil maskfile=$maskfil pbkfile=$pbkfil cleanup=no clobber=yes verbose=2
  fluximage infile=${evtfil}"[sky=region(/data/darmok/CasA/casa_box.reg)]" outroot=f016/ bands="csc,broad" xygrid="" binsize=16 asolfile=@${asollis} badpixfile=$bpixfil maskfile=$maskfil pbkfile=$pbkfil cleanup=no clobber=yes verbose=2
  fluximage infile=${evtfil}"[sky=region(/data/darmok/CasA/casa_box.reg)]" outroot=f032/ bands="csc,broad" xygrid="" binsize=32 asolfile=@${asollis} badpixfile=$bpixfil maskfile=$maskfil pbkfile=$pbkfil cleanup=no clobber=yes verbose=2
  fluximage infile=${evtfil}"[sky=region(/data/darmok/CasA/casa_box.reg)]" outroot=f064/ bands="csc,broad" xygrid="" binsize=64 asolfile=@${asollis} badpixfile=$bpixfil maskfile=$maskfil pbkfile=$pbkfil cleanup=no clobber=yes verbose=2

end
