from subprocess import call
import sys,os,time

prefix = str(sys.argv[1])
k0 = int(sys.argv[2])
k1 = int(sys.argv[3])
dk = str(sys.argv[4])
k = k0

if dk == 'None' or dk =='none':
    binFile = prefix+'.bin'
    datFile = prefix+'.dat'
    sys.stdout.write('\r {}'.format(binFile))
    sys.stdout.flush()
    call('~/bin/PolyFTS/tools/FieldBinToAscii.py -i {}  -o {}'.format(binFile,datFile),shell=True)
    while not os.path.exists(datFile):
        time.sleep(1)
    call('~/bin/PolyFTS/tools/plot/PolyFTS_to_VTK.py {}'.format(datFile),shell=True)
else:  
    dk = int(dk)     
    while k <= k1:
        binFile = prefix+'{}.bin'.format(k) 
        datFile = prefix+'{}.dat'.format(k)
        sys.stdout.write('\r {}'.format(binFile))
        sys.stdout.flush()
        call('~/bin/PolyFTS/tools/FieldBinToAscii.py -i {}  -o {}'.format(binFile,datFile),shell=True)
        while not os.path.exists(datFile):
            time.sleep(1)
        call('~/bin/PolyFTS/tools/plot/PolyFTS_to_VTK.py {}'.format(datFile),shell=True)
        k += dk

