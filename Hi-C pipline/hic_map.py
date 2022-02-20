#-*- coding: UTF-8 -*-
from sys import argv
import pandas as pd

script, suffix, name = argv

data = dict()

i = 0
for line in open(suffix+'/'+name):
  l = line.strip().split()
  k = tuple(sorted((l[2], l[6])))
  data[k] = data.get(k, 0) + 1
  i += 1
  if i % 1000000 == 0:
    print(i, flush=True)

f = open(suffix +"/" + "hic_map.txt", "w")
print("v1", "v2", "hic_w", sep="\t", file=f)
for k, v in data.items():
  print(k[0], k[1], v // 2, sep="\t", file=f)
f.close()
