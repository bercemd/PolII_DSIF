#!/usr/bin/env python

import os,sys

xlist = [0,15,30,90]
ylist = [0,15,30,90]
zlist = [0,15,30,90]
inputfile = open("5oik.seg.pdb","r")
lines = inputfile.readlines()
inputfile.close()

outputfile = open("orient.dat","w")

xcoor = []
ycoor = []
zcoor = []
for line in lines:
  if "ATOM" in line:
    xcoor.append(float(line.split()[6]))
    ycoor.append(float(line.split()[7]))
    zcoor.append(float(line.split()[8]))
box = []
box.append(max(xcoor)-min(xcoor))
box.append(max(ycoor)-min(ycoor))
box.append(max(zcoor)-min(zcoor))

boxsize = max(box)
print("boxsize",boxsize)
rotx = 0.0
roty = 0.0
rotz = 0.0
for i in range(len(xlist)):
  for j in range(len(ylist)):
    for k in range(len(zlist)):
      os.system("convpdb.pl 5oik.seg.pdb -rotatex %s | convpdb.pl -rotatey %s | convpdb.pl -rotatez %s | convpdb.pl -center > temp.pdb"%(xlist[i],ylist[j],zlist[k]))
      inputfile = open("temp.pdb","r")
      lines = inputfile.readlines()
      inputfile.close()

      xcoor = []
      ycoor = []
      zcoor = []
      for line in lines:
        if "ATOM" in line:
          xcoor.append(float(line.split()[6]))
          ycoor.append(float(line.split()[7]))
          zcoor.append(float(line.split()[8]))
      boxr = []
      boxr.append(max(xcoor)-min(xcoor))
      boxr.append(max(ycoor)-min(ycoor))
      boxr.append(max(zcoor)-min(zcoor))
      boxsizer = max(boxr)

      if boxsizer<boxsize:
        boxsize = boxsizer
        rotx = xlist[i]
        roty = ylist[j]
        rotz = zlist[k]
      print(xlist[i],ylist[j],zlist[k],boxsizer,boxsize,file=outputfile)
os.system("convpdb.pl 5oik.seg.pdb -rotatex %s | convpdb.pl -rotatey %s | convpdb.pl -rotatez %s | convpdb.pl -center > output.pdb"%(rotx,roty,rotz))
outputfile.close()
