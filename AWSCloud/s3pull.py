import boto3
import os
import sys

s3 = boto3.resource('s3')
bucket = s3.Bucket('ac-follettelab')
try:
    dirname = sys.argv[1]+'/'
except:
    dirname = ''
if(dirname  == 'options/'):
    dirname = ''

dirs = []
if(dirname == ''):
    print('options:' )
else:
    print(dirname)

for obj in bucket.objects.filter(Prefix = 'follette-lab/cloud/input/' +dirname):
    filename = obj.key.rsplit('/')

    if (dirname == ''):
        if(filename[-2] not in dirs):
            print(filename[-2])
            dirs.append(filename[-2])
    else:
        print(filename[-1])
        filename = filename[-2] + '/' + filename[-1]
        if not os.path.exists('/data/' + obj.key.rsplit('/')[-2]):
            os.makedirs('/data/' + obj.key.rsplit('/')[-2])
        s3.meta.client.download_file('ac-follettelab', obj.key, '/data/' + filename)
