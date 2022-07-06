from django.http import JsonResponse


def browse(request):
    data = {"hello": "world"}
    return JsonResponse(data)
