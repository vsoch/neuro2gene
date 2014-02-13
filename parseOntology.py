#!/usr/bin/python

import re
filey = open('12021_2013_9211_MOESM2_ESM.owl','r')
data = filey.read()


A = ">"
B = "</rdfs:label>"

matches = re.findall(re.escape(A)+"(.*?)"+re.escape(B),data)

for m in matches:
  print m
