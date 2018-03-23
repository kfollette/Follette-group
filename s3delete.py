import boto3

s3 = boto3.resource('s3')
bucket = s3.Bucket('ac-follettelab')

objects_to_delete = []
for obj in bucket.objects.filter(Prefix='follette-lab/cloud/input/' +dirname + '/'):
    objects_to_delete.append({'Key': obj.key})

bucket.delete_objects(
    Delete={
        'Objects': objects_to_delete
    }
)
