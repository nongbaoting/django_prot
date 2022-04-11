broker_url ='redis://localhost:6379/1'
# broker_url ='redis://localhost:6379/0'
# result_backend = 'redis://localhost:6379/0'
result_backend = 'django-db'
# task_serializer = 'json'
# result_serializer = 'json'
# accept_content = ['json']
enable_utc = True 
TIMEZONE='Asia/Shanghai'

# for multiprocessing, add to path
# export PYTHONOPTIMIZE=1