#!/usr/bin/csh
echo -------------------------------------------------------------------------
echo  Resubmit another run script 
echo -------------------------------------------------------------------------

set t1 = 256.25.06.15.100
set t2 = 256.25.06.15.100a
foreach file (`ls *{$t1}[.h]*`)
  echo " "
  set dmy = "`ls $file | sed -e 's/$t1/$t2/g' `"
  mv $file $dmy
end

