from os import listdir
from os.path import isfile, join
import sys

f_search_path = r'C:\Users\Wright\source\WrightTools\docs'
files = [f for f in listdir(f_search_path) if isfile(join(f_search_path, f))]

reqs = set()

for file in files:
    with open(join(f_search_path, file), 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip('\n')
            if len(line) > 0:
                l = line.split()
                if (l[0] == 'import') or (l[0] == 'from') and (l[1] != '..'):
                    reqs.add(l[1].split('.')[0])
            else:
                continue
                
print(reqs)
print([req for req in reqs if req not in sys.stdlib_module_names])



