import re
import subprocess
from collections import defaultdict
from itertools import chain
import json

def get_comparator(req):
    if req and not req[0].isalpha():
        return re.split('(\d+)',req)[0]
    else: 
        return None

def compare(a,b): 
    a_vals = [int(val) for val in a.split('.') if val.isnumeric()]
    b_vals = [int(val) for val in b.split('.') if val.isnumeric()]
    for ind in range(min([len(a_vals), len(b_vals)])):
        if a_vals[ind]==b_vals[ind]:
            continue
        elif a_vals[ind]>b_vals[ind]:
            return '>'
        elif a_vals[ind]<b_vals[ind]:
            return '<'
        else:
            print('ERROR')
    return '=='

def get_versions(lines):
    ddict = defaultdict(set)
    for line in lines.splitlines()[2:]:
        pkg, version = re.match(r"^(\S+)\s+(\S+)", line, re.MULTILINE).groups()
        ddict[pkg].add(version)
    return ddict


def semantic_cmp(version_string):
    def mysplit(string):
        version_substrs = lambda x: re.findall(r"([A-z]+|\d+)", x)
        return list(chain(map(version_substrs, string.split("."))))

    def str_ord(string):
        num = 0
        for char in string:
            num *= 255
            num += ord(char)
        return num

    def try_int(version_str):
        try:
            return int(version_str)
        except ValueError:
            return str_ord(version_str)

    mss = list(chain(*mysplit(version_string)))
    return tuple(map(try_int, mss))


def main():
    sp_i = subprocess.run(["conda", "list"], stdout=subprocess.PIPE)
    sp_v = subprocess.run(["conda", "search", "--outdated"], stdout=subprocess.PIPE)
    installed = get_versions(sp_i.stdout.decode("utf-8"))
    available = get_versions(sp_v.stdout.decode("utf-8"))
    d = dict()
    for pkg, inst_vs in installed.items():
        avail_vs = sorted(list(available[pkg]), key=semantic_cmp)
        if not avail_vs:
            continue
        current = min(inst_vs, key=semantic_cmp)
        newest = avail_vs[-1]
        if avail_vs and current != newest:
            if semantic_cmp(current) < semantic_cmp(newest):
                d[pkg]={'installed':current, 'newest':newest}
                #print(f"{pkg}:\n\tInstalled:\t{current}\n\tNewest:\t\t{newest}")
                new = [x for x in avail_vs if semantic_cmp(current) < semantic_cmp(x)]
                d[pkg]['newer'] = new
                #print("\tNewer:\t\t" + f", ".join(new))
    return d
    

if __name__ == "__main__":
    outdated = main()
    deps = subprocess.run(['pipdeptree', '--json'], stdout=subprocess.PIPE).stdout.decode('utf-8')
    deps_ls = json.loads(deps)
    '''for i in deps_ls:
        print(i['package'], '\n\n')'''
    for package in deps_ls:
        pkgname = package['package']['key']
        print('\n\n', package, '\n\n', pkgname)
        for dep in package['dependencies']:
            if dep['package_name'] in outdated: 
                current = dep['installed_version']
                reqd = dep['required_version']
                newest = outdated[dep['package_name']]['newest']
                if reqd!=None:
                    comp = get_comparator(reqd)
                    compared = compare(newest, reqd.strip(comp))
                    if comp == '>=':
                        print('\t', dep['package_name'], current, newest, reqd, compared=='>' or compared=='==')
                    elif comp == '<=':
                        print('\t', dep['package_name'], current, newest, reqd, compared=='<' or compared=='==')
                    elif comp == '==':
                        print('\t', dep['package_name'], current, newest, reqd, compared=='==')
                    elif comp == '!=':
                        print('\t', dep['package_name'], current, newest, reqd, compared!='==')