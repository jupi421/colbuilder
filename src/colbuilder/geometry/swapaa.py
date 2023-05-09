import sys
from chimera import runCommand as rc

mutation_file=str(sys.argv[1])
system_type=str(sys.argv[2])

with open(mutation_file+'.txt','r') as f:
    for l in f:
        mutation=[m for m in l.split(' ') if m!='']
        rc("open "+system_type+'/'+mutation[0])
        rc("swapaa LYS #0:"+mutation[2]+"."+mutation[3].replace('\n','').lower())
        rc("write #0 "+system_type+'/'+mutation[0])
        rc("del #0")
f.close()