function filelist,name,x

OPENR, lun1, name, /GET_LUN
list=STRARR(1000000)
a=""
WHILE (NOT EOF(lun1)) DO BEGIN
      READF, lun1, a
      list[x]=a
      x++
ENDWHILE
FREE_LUN, lun1
list=list[0:x-1]
RETURN, list

END
