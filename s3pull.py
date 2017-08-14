import boto3 
s3 = boto3.resource('s3')
bucket = s3.Bucket('ac-follettelab')
for obj in bucket.objects.all():
    filename = obj.key.rsplit(‘/‘)[-1]
    print(filename)
    s3.meta.client.download_file('ac-follettelab', obj.key, '/data/' + filename)
