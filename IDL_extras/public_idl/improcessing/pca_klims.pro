;+
; NAME: 
;   pca_klims
;
; DESCRIPTION:
;   Calculates the karhunen-loeve modes for a set of images.
;
; INPUTS:
;   err     : the covariance matrix of ims
;   rims    : a mean subtracted cube of images, in row-image format: [dim1*dim2, no_images]
;   nmodes  : the number of k-l modes to calculate
;
; INPUT KEYWORDS:
;   none
;
; KEYWORDS:
;   silent : if set, no messages are printed
;
; OUTPUTS:
;   klims   : the k-l images
;
; OUTPUT_KEYWORDS:
;    none
;
; MODIFICATION HISTORY:
;  Written 2013/01/26 by Jared R. Males (jrmales@email.arizona.edu)
;  2013/02/09 masking added by J.R.M.
;
; BUGS/WISH LIST:
;  This doesn't take full advantage of threading: neither la_eigenql nor eigenql use threadpool.
;
;-
pro pca_klims, klims, err, rims, nmodes, silent=silent

sz = size(rims)

dim1 = sz[1]
nims = sz[2]

;if(n_elements(nmodes) ne 1) then 
nmodes = nims

notsilent = 1
if(keyword_set(silent)) then notsilent = 0

if(nmodes gt nims) then begin
   if(notsilent) then begin
      print, 'pca_klims: not enough R images for ' + strcompress(string(nmodes),/rem) + ' modes, using ' + strcompress(string(nims),/rem)
   endif
   nmodes = nims
endif


;Only allocate if not passed in.  Assumes allocate was done right.
;if(n_elements(klims lt 2)) then klims = fltarr(dim1, nmodes)


if(notsilent) then begin
   print, 'pca_klims: Calculating eigenvectors . . .'
endif

;The fast way.  This returns them in ascending order, so we only want the top nmodes
lambda = la_eigenql(err,eigenvectors=eigenvectors);, range=[nims-nmodes, nims-1])

;print, lambda[-1]
if(notsilent) then begin
   print, 'pca_klims: Calculating K-L images with ' + strcompress(string(nmodes),/rem) + ' modes . . .'
endif

;Calculate the KL images:
norm = transpose(cmreplicate(1./sqrt(lambda), (size(eigenvectors))[1]))
klims = matrix_multiply(rims,norm*eigenvectors)


;la_eigenql returns eigenvs in ascending order, so reverse
klims = reverse(klims, 2, /overwrite)
   
if(notsilent) then statusline, /clear

end

