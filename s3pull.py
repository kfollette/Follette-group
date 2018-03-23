import boto3
import os

s3 = boto3.resource('s3')
bucket = s3.Bucket('ac-follettelab')
dirname = sys.argv[1]

for obj in bucket.objects.filter(Prefix = 'follette-lab/cloud/input/' +dirname + '/'):
    filename = obj.key.rsplit('/')
    filename = filename[-2] + '/' + filename[-1]
    if not os.path.exists('/data/' + obj.key.rsplit('/')[-2]):
        os.makedirs('/data/' + obj.key.rsplit('/')[-2])
    print(filename)
    s3.meta.client.download_file('ac-follettelab', obj.key, '/data/' + filename)
