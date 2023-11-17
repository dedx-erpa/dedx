from pfac import fac
import os,sys
import mprun

if __name__ == '__main__':
    np = 10
    iz0 = 0
    a0 = ''
    nz = len(mprun.zs)
    na = len(sys.argv)
    aa = 2
    ml = 2
    if na > 1:
        if sys.argv[1].isdigit():
            iz0 = int(sys.argv[1])
        else:
            a0 = sys.argv[1]
        if (len(sys.argv)>2):
            nz = int(sys.argv[2])
            if nz < 0:
                nz = len(mprun.zs)
            if (len(sys.argv)>3):
                np = int(sys.argv[3])
                if (len(sys.argv) > 4):
                    aa = int(sys.argv[4])
                    if (len(sys.argv) > 5):
                        ml = int(sys.argv[5])

    if (nz == 0):
        if iz0 > 0:
            mprun.run1zp(iz0, aa=aa, mloss=ml)
        if a0 != '':
            mprun.run1zp(a0, aa=aa, mloss=ml)
    else:
        mprun.run_azs(mprun.ss[iz0:iz0+nz], np=np, aa=aa, mloss=ml)
    
