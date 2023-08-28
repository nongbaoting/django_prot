broker_url ='redis://localhost:6379/1'
# broker_url ='redis://172.22.148.191:6379/0' # 发布
# result_backend = 'redis://localhost:6379/0' 
result_backend = 'django-db'
# task_serializer = 'json'
# result_serializer = 'json'
# accept_content = ['json']
enable_utc = True 
TIMEZONE='Asia/Shanghai'  

# for multiprocessing, add to path
# export PYTHONOPTIMIZE=1

from kombu import Exchange,Queue
task_default_queue = 'default'
task_queues = (
Queue("default",Exchange("default"),routing_key="protein.tasks.#"),
Queue("GPU",Exchange("GPU"),routing_key="protein.GPUtasks.#"),

)
# 路由
task_routes = {
'protein.GPUtasks.AlphaFold2':{"queue":"GPUtasks","routing_key":"protein.GPUtasks"},

}