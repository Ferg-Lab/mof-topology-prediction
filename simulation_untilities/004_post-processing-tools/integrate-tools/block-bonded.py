import numpy as np
from numpy import log2, zeros, mean, var, sum, loadtxt, arange, \
                  array, cumsum, dot, transpose, diagonal, floor
from numpy.linalg import inv
import sys
import os

def block(x):
    # preliminaries
    d = log2(len(x))
    if (d - floor(d) != 0):
       # print("Warning: Data size = %g, is not a power of 2." % floor(2**d))
       # print("Truncating data to %g." % 2**floor(d) )
        x = x[:2**int(floor(d))] ## length of 1025
    d = int(floor(d))
    n = 2**d
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) ) ## retreives numbers from x[0] to x[1023] * x[1] to x[1024]
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2]) ## start from 0, every 2nd elelment + start from 1, every 2nd element sum and average
        
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1] ## [::-1] reverses the array

    # we need a list of magic numbers
    q =array([6.634897,  9.210340,  11.344867, 13.276704, 15.086272, 
              16.811894, 18.475307, 20.090235, 21.665994, 23.209251,
              24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 
              31.999927, 33.408664, 34.805306, 36.190869, 37.566235,
              38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 
              45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    stdv = (s[k]/2**(d-k))**0.5
    return stdv

## calculate mean and estimate 
if len(sys.argv) < 3:
    print("usage: python block.py lambda name_of_output dirname")
    sys.exit()

lam = float(sys.argv[1])
file_name = str(sys.argv[2])
outdir = str(sys.argv[3])

os.system("grep -v '#' "+file_name+" > tmpfile && mv tmpfile "+file_name)
x = loadtxt(file_name)
time = x[:,0]
dudlam = x[:,1]

f=open('../'+outdir+'/0_bonded/averaged.dat','a')
stdv = block(dudlam)
ave = mean(dudlam)
f.write("%g %g %g \n" % (lam,ave,stdv))
f.close()

