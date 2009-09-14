#!/usr/bin/csh
echo -------------------------------------------------------------------------
echo  Resubmit another run script 
echo -------------------------------------------------------------------------

  @ N = `cat resubmit`
  if ( $N > 0 ) then
    echo "Note: resubmitting run script "
    @ N--
    echo $N >! resubmit
    bsub < run.lsf
#    cp NAMELIST1 NAMELIST
  endif

