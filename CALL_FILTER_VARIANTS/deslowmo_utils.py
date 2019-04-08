import numpy as np
import sys
import math
from platform import python_version

def get_sample(record, name):
    for samp in record.samples:
        if samp.sample == name:
            return samp
    return None

def check_record(record, c, m, f):
    mother = get_sample(record, m)
    father = get_sample(record, f)
    child = get_sample(record, c)
    return (record.is_snp and
            mother["GT"]==father["GT"] and not
            mother.is_het and 
            child.is_het)
'''
Returns True if there is a heterozygote in the record for one of the give individuals
'''
def check_hets(record, indv_names):
    indv_records = [get_sample(record, name) for name in indv_names]
    hets = [record.is_het for record in indv_records]
    return True in hets

''' *** index ***
Returns an index number based on the genotypes of the 
child: i
mother: j
father: k
based on the fact that these each have three possibilities (0,1,2)
'''
def index(i, j, k):
    return i + j*3 + k*9

''' *** trio_index ***
Returns an index number for a set of parents and their children
c_i: list of child genotypes
j: maternal genotype
k: paternal genotype
'''
def trio_index(c_i, j, k):
    result = 0
    try:
        for ind, i in enumerate(c_i):
            result += i*math.pow(3, ind)
        result += j*math.pow(3, len(c_i))
        result += k*math.pow(3, len(c_i)+1)
        result = int(result)
    except TypeError:
        print( 'child genotypes', c_i, 'are not iterable' )
    return result

''' *** ex_trio_index ***
Returns and index for an extended set of parents and offspring
c: list of child genotypes
m: list of maternal genotypes
f: list of paternal genotypes
'''
def ex_trio_index(c, m, f):
    result = 0
    try:
        _ = (x for x in enumerate(c))
        _ = (x for x in enumerate(m))
        _ = (x for x in enumerate(f))
    except TypeError:
        print( 'one of the genotype lists is not a list at all')
        return None
    for ind, i in enumerate(c):
        result += i * math.pow(3, ind)
    for ind, i in enumerate(m):
        result += i * math.pow(3, len(c) + ind)
    for ind, i in enumerate(f):
        result += i * math.pow(3, len(c) + len(m) + ind)
    return int(result)

''' *** gen_index ***
Returns an index for a sequence of genotypes
g: list of genotypes
'''
def gen_index(g):
    result = 0
    if sum([g_i > 2 for g_i in g]) > 0:
        sys.exit()
    for ind, gg in enumerate(g):
        result += gg * math.pow(3, ind)
    return int(result)

''' *** trans_matrix ***

Returns a vector the log10 probability of the child's genotype given the parental genotypes
Combinations can be retrieved using the index function
rate of new mutations per generation: mu
'''
def trans_matrix(mu):
    result = np.array([0.0]*27)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if j==0 and k==0 and i==0 :
                    result[index(i,j,k)] = (1-4*mu)*1 + (4*mu)*.5
                elif j==0 and k==1 and i==0:
                    result[index(i,j,k)] = (1-4*mu)*.5 + (4*mu)*( .5*.25 + .25*0 +.25*(.33*1 + .66*.5) )
                elif j==0 and k==0 and i==1 : 
                    result[index(i,j,k)] = (1-4*mu)*0 + (4*mu)*(.33*.5 + .66*0)
                elif j==0 and k==1 and i==1 : 
                    result[index(i,j,k)] = (1-4*mu)*.5 + (4*mu)*(.5*(.33*.5 + .66*.5*.5) + .25*(.33*1 + .66*.5) + .25*0)
                elif j==0 and k==2 and i==0 : 
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.5*0 + .5*(.33*.5 + .66*0))
                elif j==0 and k==2 and i==1 :
                    result[index(i,j,k)] =  (1-4*mu)*1 + (4*mu)*(.5*.5 + .5*.5)
                elif j==0 and k==2 and i==2 :
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.5*(.33*.5 + .66*0) + .5*0)
                elif j==0 and k==0 and i==2 :
                    result[index(i,j,k)] =  sys.float_info.min
                elif j==0 and k==1 and i==2 :
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.5*(.33*.5*.5 + .66*0) + .5*0)
                elif j==1 and k==0 and i==0 :
                    result[index(i,j,k)] =  (1-4*mu)*.5 + (4*mu)*(.25*(.33*1+.66*.5) + .25*0 + .5*.25)
                elif j==1 and k==1 and i==0 :
                    result[index(i,j,k)] =  (1-4*mu)*.25 + (4*mu)*(.5*0 + .5*(.33*.5 + .66*.25))
                elif j==1 and k==0 and i==1 :
                    result[index(i,j,k)] =  (1-4*mu)*.5 + (4*mu)*(.5*(.33*.5 + .66*.5*.5) + .25*(.33*1 + .66*.5) + .25*0)
                elif j==1 and k==1 and i==1 :
                    result[index(i,j,k)] =  (1-4*mu)*.5 + (4*mu)*(.5*(.33*.5+.66*.5*.5) + .5*(.33*.5+.66*.5*.5))
                elif j==1 and k==2 and i==0 :
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.5*0 + .5*(.33*.5*.5 + .66*0))
                elif j==1 and k==2 and i==1 :
                    result[index(i,j,k)] =  (1-4*mu)*.5 + (4*mu)*(.5*(.33*.5 + .66*.5*.5) + .25*(.33*1 + .66*.5) + .25*0)
                elif j==1 and k==2 and i==2 :
                    result[index(i,j,k)] =  (1-4*mu)*.5 + (4*mu)*(.25*(.33*1+.66*.5) + .25*0 + .5*.25)
                elif j==1 and k==0 and i==2 :
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.5*(.33*.5*.5 + .66*0) + .5*0)
                elif j==1 and k==1 and i==2 :
                    result[index(i,j,k)] =  (1-4*mu)*.25 + (4*mu)*(.5*0 + .5*(.33*.5 + .66*.25))
                elif j==2 and k==0 and i==0 :
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.5*(.33*.5 + .66*0) + .5*0)
                elif j==2 and k==1 and i==0 :
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.5*(.33*.5*.5 + .66*0) + .5*0)
                elif j==2 and k==0 and i==1 :
                    result[index(i,j,k)] =  (1-4*mu)*1 + (4*mu)*(.5*.5 + .5*.5)
                elif j==2 and k==1 and i==1 :
                    result[index(i,j,k)] =  (1-4*mu)*.5 + (4*mu)*(.5*(.33*.5 + .66*.5*.5) + .25*(.33*1 + .66*.5) + .25*0)
                elif j==2 and k==2 and i==0 :
                    result[index(i,j,k)] =  sys.float_info.min
                elif j==2 and k==2 and i==1 :
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.33*.5 + .66*0)
                elif j==2 and k==2 and i==2 :
                    result[index(i,j,k)] =  (1-4*mu)*1 + (4*mu)*.5
                elif j==2 and k==0 and i==2 :
                    result[index(i,j,k)] =  (1-4*mu)*0 + (4*mu)*(.5*0 + .5*(.33*.5 + .66*0))
                elif j==2 and k==1 and i==2 :
                    result[index(i,j,k)] =  (1-4*mu)*.5 + (4*mu)*( .5*.25 + .25*0 +.25*(.33*1 + .66*.5) )
    return np.log10(result)
