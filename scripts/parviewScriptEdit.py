#!/usr/bin/env python3

#Alters Pic taking scripts with rescaling factors and date


import re
import argparse
import fileinput
from sys import argv

# Prep for image taking 

clv3ARescale="RescaleTransferFunction"+ argv[2]

clv3BRescale="RescaleTransferFunction"+ argv[3]

clv3PRescale="RescaleTransferFunction"+ argv[4]

dimerRescale="RescaleTransferFunction"+ argv[5]

monomerRescale="RescaleTransferFunction"+ argv[6]

wusCRescale="RescaleTransferFunction"+ argv[7]

wusNRescale="RescaleTransferFunction"+ argv[8]

wusRRescale="RescaleTransferFunction"+ argv[9]

inputfile= 'weitaoextended7.py'

# Read in the file
with open( inputfile, 'r') as file :
  filedata = file.read()


# Replace the target string
filedata = filedata.replace('RescaleTransferFunction(0.0, 1.0)', clv3ARescale)

# Replace the target string
filedata = filedata.replace('RescaleTransferFunction(0.0, 2.0)', clv3BRescale)

# Replace the target string
filedata = filedata.replace('RescaleTransferFunction(0.0, 3.0)', clv3PRescale)

# Replace the target string
filedata = filedata.replace('RescaleTransferFunction(0.0, 4.0)', dimerRescale)

# Replace the target string
filedata = filedata.replace('RescaleTransferFunction(0.0, 5.0)', monomerRescale)

# Replace the target string
filedata = filedata.replace('RescaleTransferFunction(0.0, 6.0)', wusCRescale)

# Replace the target string
filedata = filedata.replace('RescaleTransferFunction(0.0, 7.0)', wusNRescale)

# Replace the target string
filedata = filedata.replace('RescaleTransferFunction(0.0, 8.0)', wusRRescale)

# Replace the target string
filedata = filedata.replace('021121DcTrWt', argv[1])
filedata = filedata.replace('021121FcAcWt', argv[1])

outputfile= 'weitaoextended7Temp.py'

# Write the file out again
with open(outputfile, 'w+') as file:
  file.write(filedata)







