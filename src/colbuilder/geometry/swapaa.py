import sys
from chimera import runCommand as rc

file=str(sys.argv[1])
system_type=str(sys.argv[2])

with open(file+'.txt','r') as f:
    for l in f:
        replace=[m for m in l.split(' ') if m!='']
        rc("open "+system_type+'/'+replace[0])
        rc("swapaa LYS #0:"+replace[2]+"."+replace[3].replace('\n','').lower())
        rc("write #0 "+system_type+'/'+replace[0])
        rc("del #0")
f.close()