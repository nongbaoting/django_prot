import re
import time
import json
from protein.models import *
from django.middleware import csrf
from collections import defaultdict
from django.utils import timezone
from django.shortcuts import render
from django.core.files import File
from django.http import JsonResponse, HttpResponse
from django.core import serializers
import datetime

userDt = {'username': 'teamwork',
'password': 'weareteam!'}
def login(request):
    body_unicode = request.body.decode('utf-8')
    body = json.loads(body_unicode)
    user = body['username']
    passwd = body['password']
    data = {"status":-1}
    print(user, passwd)
    token = csrf.get_token(request)
    if user == userDt['username'] and passwd == userDt['password']:
        data['status'] = 200
        data['token'] = token
    return JsonResponse(data)