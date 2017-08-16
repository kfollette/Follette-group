import boto3
import os
import sys

s3 = boto3.resource('s3')
dirname = sys.argv[1]

for f in os.listdir('/data/' + dirname):
    print(f)
    s3.Object('ac-follettelab', 'follette-lab/cloud/output/' + dirname + '/' + f).upload_file('/data/' + dirname + '/' + f)
