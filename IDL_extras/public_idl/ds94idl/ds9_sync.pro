pro ds9_sync

common ds9env, $    ;this is the ds9 environment common block
       ds9obj    ;the ds9 object
       


;check if ds9 is initialized
ds9_init

ds9obj->dumpfile, fn

ds9obj->setimage, mrdfits(fn)

end


