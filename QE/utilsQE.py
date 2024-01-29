import numpy as np
import matplotlib.pyplot as plt

def read_filband(file='filband'):
    filelines = [line for line in open(file) if line.strip()]
    
    nbnd = int(filelines[0].split()[2].split(',')[0])
    nks = int(filelines[0].split()[4])
    
    kpointline = True
    for i, line in enumerate(filelines[1:]):
        if kpointline:
            try:
                kpoints=np.vstack([kpoints,np.array(line.split(), dtype='float')])
            except NameError:
                kpoints=np.array(np.array(line.split(), dtype='float'))
            
            bandval = []
            kpointline=False
        else:
            for dummy in line.split():
                bandval.append(float(dummy))
            if len(bandval) == nbnd:
                try:
                    bands=np.vstack([bands,np.array(bandval)])
                except NameError:
                    bands=np.array(bandval)
                kpointline=True
    
    return kpoints, bands


def read_bandsin(bandsinfile='bands.in'):
    filelines = [line for line in open(bandsinfile) if line.strip()]
    
    inKpoints=False
    inCellPar = False
    highsympoints = []
    pointsinline=[]
    
    for i, line in enumerate(filelines):
        if not inCellPar:
            if line.split()[0].upper() == 'CELL_PARAMETERS':
                inCellPar = True
                counter = 0
        else:
            try:
                B=np.vstack([B,np.array(line.split()[-3:],dtype=float)])
            except NameError:
                B=np.array(np.array(line.split()[-3:],dtype=float))
            counter += 1
            if counter == 3: 
                B = np.array(B)
                inCellPar=False
            
            
        
        if not inKpoints:
            if line.split()[0].upper() == 'K_POINTS':
                inKpoints = True
        else:
            if len(line.split()) == 1:
                nlines = int(line)
            else:
                try:
                    kpoints=np.vstack([kpoints,np.array(line.split()[0:3], dtype='float')])
                except NameError:
                    kpoints=np.array(np.array(line.split()[0:3], dtype='float'))
                pointsinline.append(int(line.split()[3]))
                hsp = line.split('!')[1].split('\n')[0]
                if hsp.upper() == 'GAMMA': hsp = '\Gamma'
                highsympoints.append('${}$'.format(hsp))
                
                if len(highsympoints) == nlines: inKpoints=False
    
    return kpoints, nlines, highsympoints, B, pointsinline


def readPhononbandFreq(freqgpfile): #freqgpfile=SYSTEM.freq.gp

    data = np.loadtxt(freqgpfile)

    nbands = data.shape[1] - 1
    qpoints = data[:, 0]
    bands = np.transpose(data)[1:,]

    return qpoints, bands


def readHighSymPointsPhonon(matdynfile):  #matdynfile=matdyn.in
    filelines = [line for line in open(matdynfile) if line.strip()]
    qbandform=False
    for i,fline in enumerate(filelines):
        
        if fline.split()[0]=='readtau' and fline.split()[-1]=='.true.':
            continue #TODO: INCLUDE THE CASE WHERE ATOMIC POSITIONS ARE READ

        if fline.split()[0]=='q_in_band_form' and fline.split()[-1]=='.true.':
            qbandform=True

        if fline.split()[0]=='/':
            if qbandform:
                start = i+1

            #TODO: make the case where qbandform=False
    
    npoints = int(filelines[start])
    labels = []
    positions = np.empty(npoints, dtype='int')
    counter = 0
    for i,kline in enumerate(filelines[start+1:]):
        dummylabel = str(kline.split('!')[-1].split()[0])
        if dummylabel[0].upper() == 'G':
            labels.append('$\\Gamma$')  
        else: labels.append('${}$'.format(dummylabel))
        positions[i] = counter
        counter += int(kline.split()[-2])
    return labels, positions