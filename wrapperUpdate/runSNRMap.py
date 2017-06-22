
# coding: utf-8

# In[1]:

import SNRMap as SNR
import sys
############ MODULE TESTING #################
#print2D(implant_custom_mask(45,90,150,200))
#############################################

#########################################################################################################################################################################################################
###########################################################################       RUN PROGRAM        ####################################################################################################
file = sys.argv[1]
output = sys.argv[2]
new = sys.argv[3]
if new == 'True':
    new = True
elif new == 'False':
    new = False
graph = sys.argv[4]
if graph == 'True':
    graph = True
elif graph == 'False':
    graph = False
if len(sys.argv) > 5:
    theta1 = float(sys.argv[5])
    theta2 = float(sys.argv[6])
    r1 = float(sys.argv[7])
    r2 = float(sys.argv[8])
    mask = SNR.implant_custom_mask(theta1, theta2, r1, r2)
    test = SNR.SNRMap(file, output, new, graph, mask)
else:
    test = SNR.SNRMap(file, output, new, graph)
#med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits
#########################################################################################################################################################################################################



# In[ ]:



