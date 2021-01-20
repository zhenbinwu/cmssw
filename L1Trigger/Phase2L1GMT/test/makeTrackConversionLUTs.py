import math
BITSABSCURV=12-1
BITSPT=13
maxCurv = 0.00855
ptLSB=0.025
ptLUT=[]
for i in range(1,(1<<BITSABSCURV)):
    k = (maxCurv*i)/(1<<BITSABSCURV)
    pOB=0.3*3.8*0.01/(k)
    pINT = int(pOB/ptLSB)
    if pINT<(1<<BITSPT):
        print i,k,pOB,pINT
        ptLUT.append(str(pINT))
    else:    
        ptLUT.append(str((1<<BITSPT)-1))


print("ap_uint<BITSPT> ptLUT[{address}]={{".format(address=(1<<BITSABSCURV))+','.join(ptLUT)+'};')



BITSTTTANL=16-1-4
BITSETA=12
maxTanL=8.0
etaLUT=[]
for i in range(0,(1<<(BITSTTTANL))):
    tanL = (maxTanL*i)/(1<<BITSTTTANL)
    lam =math.atan(tanL)
    theta =math.pi/2.0-lam 
    eta = -math.log(math.tan(theta/2.0))
    etaINT = int(eta*(1<<BITSETA)/math.pi)
    if abs(eta<math.pi):
        etaLUT.append(str(etaINT))
#        etaLUT.append('ap_uint<BITSETA-1>({val})'.format(val=etaINT))


print("ap_uint<BITSETA-1> etaLUT[{address}]={{".format(address=(1<<BITSTTTANL))+','.join(etaLUT)+'};')


