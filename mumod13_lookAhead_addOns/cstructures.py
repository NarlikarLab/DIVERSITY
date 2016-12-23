from ctypes import *

class DataSet(Structure):
    _fields_ = [("data", POINTER(POINTER(c_int))), ("features", POINTER(c_int)), ("lookAhead", POINTER(POINTER(c_int))), ("featureValues", c_int), ("n", c_int), ("tu", c_int)]

class TrainOut(Structure):
    _fields_ = [("labels", POINTER(c_int)), ("startPos", POINTER(c_int)), ("motifWidth", POINTER(c_int)), ("likelihood", c_double)]

def getCArray(l, c):  # Convert Python list to C ID array
    n = len(l)
    arrOfN = c*n
    num = arrOfN()
    for i in xrange(n): num[i] = l[i]
    return num

def getTrainOut(to, mode, n):
    d = {}
    d['labels'] = [to.contents.labels[i] for i in xrange(n)]
    d['startPos'] = [to.contents.startPos[i] for i in xrange(n)]
    d['motifWidth'] = [to.contents.motifWidth[i] for i in xrange(mode)]
    d['likelihood'] = to.contents.likelihood
    return d

libctest = cdll.LoadLibrary("/denalinfs/home/home/anushua/mumod13_lookAhead_addOns/libdatacal.so")
libctest.getData.restype = POINTER(DataSet)
libctest.getBackground1.restype = POINTER(POINTER(c_double))
libctest.trainData.restype = POINTER(TrainOut)
