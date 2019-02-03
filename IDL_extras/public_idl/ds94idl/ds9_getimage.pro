function ds9_getimage

common ds9env, $    ;this is the ds9 environment common block
       ds9obj    ;the ds9 object
       


ds9obj->getimage, im

return, im

end
