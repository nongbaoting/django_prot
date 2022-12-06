def cas12f(request):


def parse_info(infoFile):
    data = []
    with open(infoFile, 'r') as f:
        header = next(f).strip('\n').split("\t")
        header_len = len(header)
        for li in f:
            cell = li.split('\n').split("\t")
            # dt = {
            #     header[i]: cell[i] for i in range(header_len)
            # }
            dt = {

            }

        data.append(dt)
    return data
