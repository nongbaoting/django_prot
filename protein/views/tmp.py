# 分页

# paginator = Paginator(newQ, pageSize)
# print("currentPage, pageSize",currentPage, pageSize)
# try:
#     books = paginator.page(currentPage)
# except PageNotAnInteger:
#     books = paginator.page(1)
# except EmptyPage:
#     books = paginator.page(paginator.num_pages)
# content = serializers.serialize('json',books)

      
# qset = []
# for domain in body['domains']:
#     if domain['value'] and domain['type']=='cdd_name' and not domain["exclude"]:
#         dCount = domain['count']
#         keyword =  domain['value'].strip()
#         print("key word", keyword, "count: ", dCount)
#         re_k = re.compile(f"({keyword})", re.I)
#         for query in querySet:
#             print(query.cdd_annots)
#             matchs = len(re_k.findall(query.cdd_annots))
#             if matchs >= dCount:
#                 qset.append(query)

# print(len(qset))
# print(json.dumps(qset))