import sys 
def makingFullFilePath_linux(cwd=None, file=None): 
    if isinstance(cwd, type(None)): 
        cwd = sys.getcwd()
    if cwd[-1] != "/": cwd+='/'

    if "../" in file: 
        names  =file.strip().split("../")
        cnt = 0 
        for name in names: 
            if "" in name: 
                cnt += 1 
        dirs = cwd.split("/")
        path =""
        for i in range(len(dirs)-cnt): 
            if dirs[i]=='':
                path +='/'
            else:
                path += dirs[i]+'/'
        for name in names: 
            if name != "": 
                path += name 
    elif '/home' == file[:5]: 
        path = file 
    else: 
        if '/' == file[0]: 
            file = file[1:]
        path = cwd + file 

    print (path)
    return path 

def makeUpforShortage(lst=None, N=0, replace=None): 
    if isinstance(replace, type(None)): 
        replace = lst[-1]
    
    i = len(lst)
    while i < N: 
        lst.append(replace)
        i += 1 
    return lst 

def split_string(string=None, divider=None): 
    if isinstance(divider, type(None)): 
        return [string]
    wds = string.split(divider)
    words =[]
    for wd in wds: 
        if wd.strip() != '': 
            words.append(wd.strip())
    return words 

def series(word): 
    wds = word.split("~")
    wd1 = int(wds[0].strip())
    wd2 = int(wds[1].strip())
    srs = []
    for i in range(wd1, wd2+1): 
        srs.append(str(i))
    return srs 

def simulationCodes(vts, vnums, revs, types, snums, doe=False): 
    vts = split_string(vts, '/')
    vnums= split_string(vnums, '/')
    revs = split_string(revs, '/')
    types = split_string(types, '/')
    snums = split_string(snums, '/')

    N = len(vts)
    vnums = makeUpforShortage(vnums, N, replace='0')
    revs = makeUpforShortage(revs, N, replace='0')
    types = makeUpforShortage(types, N, replace='D101')
    snums = makeUpforShortage(snums, N, replace='0')


    sims=[]
    for vt, vnum, rev, tp, snum in zip(vts, vnums, revs, types, snums): 
        tmp =[vt]
        # print (vt, ',', vnum, ',', rev, ',', tp, ',', snum)
        vns = split_string(vnum, ",")
        rvs = split_string(rev, ',')
        tps = split_string(tp, ',')
        sns = split_string(snum, ',')
        for vn in vns: 
            if '~' in vn:  vn = series(vn)
            else: vn =[vn]
            for v in vn: 
                for rv in rvs: 
                    if '~' in rv:  rv = series(rv)
                    else: rv = [rv]
                    for r in rv: 
                        for ts in tps: 
                            for sn in sns: 
                                if '~' in sn:  sn = series(sn)
                                else: sn =[sn]
                                for s in sn : 
                                    sims.append([vt, v, r, ts, s])
    return sims 


if __name__ =="__main__": 
    cwd = '/home/users/h202/2021/'
    file = '../../aaa/bbb/ccc.txt'
    file = '/aaa/bbb.txt'
    file = 'aa/bbb.txt'
    file = '../home/txt'
    makingFullFilePath_linux(cwd, file)
