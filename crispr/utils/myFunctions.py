# def scanAndFind_pattern(mydir, mypattern):
#     wantFiles = []
#     for entry in os.scandir(mydir):
#         if entry.is_file() and mypattern.search(entry.name):
#             wantFiles.append(entry)
#         elif entry.is_dir():
#             wantFiles.extend(scanAndFind_pattern(entry.path, mypattern))
#         elif entry.is_symlink():
#             mydir = os.path.abspath(entry.path)
#             if os.path.isdir(mydir):
#                 print(mydir)
#                 wantFiles.extend(scanAndFind_pattern(
#                     mydir, mypattern))
#     return wantFiles
import os
import subprocess


def scanAndFind_pattern(mydir, mypattern):
    wantFiles = []
    for entry in os.scandir(mydir):
        if entry.is_file() and mypattern.search(entry.name):
            wantFiles.append(entry)
        elif entry.is_dir() or entry.is_symlink():
            wantFiles.extend(scanAndFind_pattern(entry.path, mypattern))
    return wantFiles


def calAlignIdent(a, b, c):
    eff_indices = [index for index, item in enumerate(c) if item in [":"]]
    align_indices = [index for index, item in enumerate(c) if item in [
        ":", "."]]

    align_len = 0
    for i in align_indices:
        if a[i] == b[i]:
            align_len += 1
    align4_len = 0

    for i in eff_indices:
        if a[i] == b[i]:
            align4_len += 1

    align_ident = round(align_len/len(align_indices), 2)

    return [align4_len, align_ident]


def run_cmd(cmd):
    run = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return run
