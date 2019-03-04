#!/usr/bin/python
import os

path = './29.4.2016_RACK1siRNA-BactinFISH'
dirs = os.listdir(path)
for dir in dirs:
    command = "python create_3d_image.py -t FISH"
    command = command+" -o "+os.path.join(path, dir).replace(" ", "\ ")
    command = command + " -i " +os.path.join(path, dir).replace(" ", "\ ")
    print(command)
    os.system(command)
