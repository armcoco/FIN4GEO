from t2listing import *
lst=t2listing('outputShorterFaults/flow.out')
xcenter=np.loadtxt('outputShorterFaults/xcenter.txt')
ycenter=np.loadtxt('outputShorterFaults/ycenter.txt')
zcenter=np.loadtxt('outputShorterFaults/zcenter.txt')
for iter in range(0, lst.num_fulltimes):
    lst.index=iter
    Temp=lst.element['T']
    Pres=lst.element['P']
    SG=lst.element['SG']
    DW=lst.element['DW']
    DG=lst.element['DG']
    np.savetxt('outputShorterFaults/Sol_'+str(iter),np.c_[xcenter,ycenter,zcenter,Temp,Pres,SG,DW,DG])

np.savetxt('outputShorterFaults/time_heat',np.c_[lst.times])


