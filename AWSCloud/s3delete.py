import boto3
import sys
import glob

s3 = boto3.resource('s3')
bucket = s3.Bucket('amherst-follette-lab')
dirname = sys.argv[1]

if dirname=='all':
	go = input('Are you sure you would like to delete everything in the s3? (y/n)')
	if go=='y' or go=='Y':
		bucket.object_versions.delete()
	else:
		print('ok. doing nothing.')

else:
	objects_to_delete = []
	for obj in bucket.objects.filter(Prefix='follette-lab/cloud/input/' +dirname + '/'):
	    objects_to_delete.append({'Key': obj.key})


	try:
	    bucket.delete_objects(
	        Delete={'Objects': objects_to_delete})
	except:
	    for obj in bucket.objects.filter(Prefix='follette-lab/cloud/input/' +dirname + '/'):
	        obj.delete()
