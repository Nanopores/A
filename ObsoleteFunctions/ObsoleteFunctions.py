def SaveToHDF5(inp_file, AnalysisResults, coefficients, outdir):
    file = str(os.path.split(inp_file['filename'])[1][:-4])
    f = h5py.File(outdir + file + '_OriginalDB.hdf5', "w")
    general = f.create_group("General")
    general.create_dataset('FileName', data=inp_file['filename'])
    general.create_dataset('Samplerate', data=inp_file['samplerate'])
    general.create_dataset('Machine', data=inp_file['type'])
    general.create_dataset('TransverseRecorded', data=inp_file['graphene'])
    general.create_dataset('TimeFileWritten', data=inp_file['TimeFileWritten'])
    general.create_dataset('TimeFileLastModified', data=inp_file['TimeFileLastModified'])
    general.create_dataset('ExperimentDuration', data=inp_file['ExperimentDuration'])

    segmentation_LP = f.create_group("LowPassSegmentation")
    for k,l in AnalysisResults.items():
        set1 = segmentation_LP.create_group(k)
        lpset1 = set1.create_group('LowPassSettings')
        for o, p in coefficients[k].items():
             lpset1.create_dataset(o, data=p)
        for m, l in AnalysisResults[k].items():
            if m is 'AllEvents':
                eventgroup = set1.create_group(m)
                for i, val in enumerate(l):
                    eventgroup.create_dataset('{:09d}'.format(i), data=val)
            elif m is 'Cusum':
                eventgroup = set1.create_group(m)
                for i1, val1 in enumerate(AnalysisResults[k]['Cusum']):
                    cusevent = eventgroup.create_group('{:09d}'.format(i1))
                    cusevent.create_dataset('NumberLevels', data=np.uint64(len(AnalysisResults[k]['Cusum'][i1]['levels'])))
                    if len(AnalysisResults[k]['Cusum'][i1]['levels']):
                        cusevent.create_dataset('up', data=AnalysisResults[k]['Cusum'][i1]['up'])
                        cusevent.create_dataset('down', data=AnalysisResults[k]['Cusum'][i1]['down'])
                        cusevent.create_dataset('both', data=AnalysisResults[k]['Cusum'][i1]['both'])
                        cusevent.create_dataset('fit', data=AnalysisResults[k]['Cusum'][i1]['fit'])
                        # 0: level number, 1: current, 2: length, 3: std
                        cusevent.create_dataset('levels_current', data=AnalysisResults[k]['Cusum'][i1]['levels'][1])
                        cusevent.create_dataset('levels_length', data=AnalysisResults[k]['Cusum'][i1]['levels'][2])
                        cusevent.create_dataset('levels_std', data=AnalysisResults[k]['Cusum'][i1]['levels'][3])
            else:
                set1.create_dataset(m, data=l)

def CorrelateTheTwoChannels(AnalysisResults, DelayLimit, Ch1 = 'i1', Ch2 = 'i2'):
    if len(AnalysisResults[Ch1]['RoughEventLocations']) is not 0:
        i1StartP = np.int64(AnalysisResults[Ch1]['StartPoints'][:])
    else:
        i1StartP = []
    if len(AnalysisResults[Ch2]['RoughEventLocations']) is not 0:
        i2StartP = np.int64(AnalysisResults[Ch2]['StartPoints'][:])
    else:
        i2StartP = []

    # Common Events, # Take Longer
    CommonEventsi1Index = np.array([], dtype=np.uint64)
    CommonEventsi2Index = np.array([], dtype=np.uint64)

    for k in i1StartP:
        val = i2StartP[(i2StartP > k - DelayLimit) & (i2StartP < k + DelayLimit)]
        if len(val) == 1:
            CommonEventsi2Index = np.append(CommonEventsi2Index, np.uint64(np.where(i2StartP == val)[0][0]))
            CommonEventsi1Index = np.append(CommonEventsi1Index, np.uint64(np.where(i1StartP == k)[0][0]))
        if len(val) > 1:
            diff = np.absolute(val-k)
            CommonEventsi2Index = np.append(CommonEventsi2Index, np.uint64(np.where(i2StartP == val[np.argmin(diff)]))[0][0])
            CommonEventsi1Index = np.append(CommonEventsi1Index, np.uint64(np.where(i1StartP == k)[0][0]))
    # Only i1
    Onlyi1Indexes = np.delete(range(len(i1StartP)), CommonEventsi1Index)
    #Onlyi1Indexes=[]
    # Only i2
    Onlyi2Indexes = np.delete(range(len(i2StartP)), CommonEventsi2Index)
    #Onlyi2Indexes=[]

    CommonIndexes={}
    CommonIndexes[Ch1]=CommonEventsi1Index
    CommonIndexes[Ch2]=CommonEventsi2Index
    OnlyIndexes={}
    OnlyIndexes[Ch1] = Onlyi1Indexes
    OnlyIndexes[Ch2] = Onlyi2Indexes
    return (CommonIndexes, OnlyIndexes)

def PlotEvent(fig1, t1, i1, t2 = [], i2 = [], fit1 = np.array([]), fit2 = np.array([]), channel = 'i1'):
    if len(t2)==0:
        ax1 = fig1.add_subplot(111)
        ax1.plot(t1, i1*1e9, 'b')
        if len(fit1) is not 0:
            ax1.plot(t1, fit1*1e9, 'y')
        ax1.set_ylabel(channel + ' Current [nA]')
        ax1.set_xlabel(channel + ' Time [s]')
        ax1.ticklabel_format(useOffset=False)
        ax1.ticklabel_format(useOffset=False)
        return ax1
    else:
        ax1 = fig1.add_subplot(211)
        ax2 = fig1.add_subplot(212, sharex=ax1)
        ax1.plot(t1, i1*1e9, 'b')
        if len(fit1) is not 0:
            ax1.plot(t1, fit1*1e9, 'y')
        ax2.plot(t2, i2*1e9, 'r')
        if len(fit2) is not 0:
            ax2.plot(t2, fit2*1e9, 'y')
        ax1.set_ylabel('Ionic Current [nA]')
        #ax1.set_xticklabels([])
        ax2.set_ylabel('Transverse Current [nA]')
        ax2.set_xlabel('Time [s]')
        ax2.ticklabel_format(useOffset=False)
        ax2.ticklabel_format(useOffset=False)
        ax1.ticklabel_format(useOffset=False)
        ax1.ticklabel_format(useOffset=False)
        return ax1, ax2

def SaveAllPlots(CommonIndexes, OnlyIndexes, AnalysisResults, directory, out, buffer, withFit = 1):
    if len(CommonIndexes['i1']) is not 0:
        # Plot All Common Events
        pp = PdfPages(directory + '_SavedEventsCommon.pdf')
        ind1 = np.uint64(CommonIndexes['i1'])
        ind2 = np.uint64(CommonIndexes['i2'])

        t = np.arange(0, len(out['i1']))
        t = t / out['samplerate'] * 1e3
        count=1
        for eventnumber in range(len(ind1)):
            parttoplot = np.arange(AnalysisResults['i1']['StartPoints'][ind1[eventnumber]] - buffer,
                                   AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] + buffer, 1, dtype=np.uint64)
            parttoplot2 = np.arange(AnalysisResults['i2']['StartPoints'][ind2[eventnumber]] - buffer,
                                    AnalysisResults['i2']['EndPoints'][ind2[eventnumber]] + buffer, 1, dtype=np.uint64)

            fit1 = np.concatenate([np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]],
                                   np.ones(AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] - AnalysisResults['i1']['StartPoints'][
                                       ind1[eventnumber]]) * (
                                       AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]] - AnalysisResults['i1']['DeltaI'][ind1[eventnumber]]),
                                   np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]]])

            fit2 = np.concatenate([np.ones(buffer) * AnalysisResults['i2']['LocalBaseline'][ind2[eventnumber]],
                                   np.ones(AnalysisResults['i2']['EndPoints'][ind2[eventnumber]] - AnalysisResults['i2']['StartPoints'][
                                       ind2[eventnumber]]) * (
                                       AnalysisResults['i2']['LocalBaseline'][ind2[eventnumber]] - AnalysisResults['i2']['DeltaI'][ind2[eventnumber]]),
                                   np.ones(buffer) * AnalysisResults['i2']['LocalBaseline'][ind2[eventnumber]]])
            if withFit:
                fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot2], out['i2'][parttoplot2],
                                   fit1=fit1, fit2=fit2)
            else:
                fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot2], out['i2'][parttoplot2])

            if not np.mod(eventnumber+1,200):
                pp.close()
                pp = PdfPages(directory + '_SavedEventsCommon_' + str(count) + '.pdf')
                count+=1
            pp.savefig(fig)
            print('{} out of {} saved!'.format(str(eventnumber), str(len(ind1))))
            print('Length i1: {}, Fit i1: {}'.format(len(out['i1'][parttoplot]), len(fit1)))
            print('Length i2: {}, Fit i2: {}'.format(len(out['i2'][parttoplot2]), len(fit2)))
            fig.clear()
            plt.close(fig)
        pp.close()

    if len(OnlyIndexes['i1']) is not 0:
        # Plot All i1
        pp = PdfPages(directory + '_SavedEventsOnlyi1.pdf')
        ind1 = np.uint64(OnlyIndexes['i1'])

        t = np.arange(0, len(out['i1']))
        t = t / out['samplerate'] * 1e3
        count=1
        for eventnumber in range(len(ind1)):
            parttoplot = np.arange(AnalysisResults['i1']['StartPoints'][ind1[eventnumber]] - buffer,
                                   AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] + buffer, 1, dtype=np.uint64)

            fit1 = np.concatenate([np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]],
                                   np.ones(AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] - AnalysisResults['i1']['StartPoints'][
                                       ind1[eventnumber]]) * (
                                       AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]] - AnalysisResults['i1']['DeltaI'][ind1[eventnumber]]),
                                   np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]]])

            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot], out['i2'][parttoplot], fit1=fit1)
            if not np.mod(eventnumber+1,200):
                pp.close()
                pp = PdfPages(directory + '_SavedEventsCommon_' + str(count) + '.pdf')
                count+=1
            pp.savefig(fig)
            print('{} out of {} saved!'.format(str(eventnumber), str(len(ind1))))
            fig.clear()
            plt.close(fig)
        pp.close()

    if len(OnlyIndexes['i2']) is not 0:
        # Plot All i2
        pp = PdfPages(directory + '_SavedEventsOnlyi2.pdf')
        ind1 = np.uint64(OnlyIndexes['i2'])

        t = np.arange(0, len(out['i2']))
        t = t / out['samplerate'] * 1e3
        count=1
        for eventnumber in range(len(ind1)):
            parttoplot = np.arange(AnalysisResults['i2']['StartPoints'][ind1[eventnumber]] - buffer,
                                   AnalysisResults['i2']['EndPoints'][ind1[eventnumber]] + buffer, 1, dtype=np.uint64)

            fit1 = np.concatenate([np.ones(buffer) * AnalysisResults['i2']['LocalBaseline'][ind1[eventnumber]],
                                   np.ones(AnalysisResults['i2']['EndPoints'][ind1[eventnumber]] - AnalysisResults['i2']['StartPoints'][
                                       ind1[eventnumber]]) * (
                                       AnalysisResults['i2']['LocalBaseline'][ind1[eventnumber]] - AnalysisResults['i2']['DeltaI'][ind1[eventnumber]]),
                                   np.ones(buffer) * AnalysisResults['i2']['LocalBaseline'][ind1[eventnumber]]])

            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot], out['i2'][parttoplot], fit2=fit1)
            if not np.mod(eventnumber+1,200):
                pp.close()
                pp = PdfPages(directory + '_SavedEventsCommon_' + str(count) + '.pdf')
                count+=1
            pp.savefig(fig)
            print('{} out of {} saved!'.format(str(eventnumber), str(len(ind1))))
            fig.clear()
            plt.close(fig)
        pp.close()

    # Derivative
    if len(CommonIndexes['i1']) is not 0:
        # Plot All i1
        pp = PdfPages(directory + '_i1vsderivi2.pdf')
        ind1 = np.uint64(CommonIndexes['i1'])
        ind2 = np.uint64(CommonIndexes['i2'])

        t = np.arange(0, len(out['i1']))
        t = t / out['samplerate'] * 1e3
        count=1
        for eventnumber in range(len(ind1)):
            parttoplot = np.arange(AnalysisResults['i1']['StartPoints'][ind1[eventnumber]] - buffer,
                                   AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] + buffer, 1, dtype=np.uint64)
            parttoplot2 = np.arange(AnalysisResults['i2']['StartPoints'][ind2[eventnumber]] - buffer,
                                    AnalysisResults['i2']['EndPoints'][ind2[eventnumber]] + buffer, 1, dtype=np.uint64)

            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot2][:-1],
                               np.diff(out['i2'][parttoplot2]))

            if not np.mod(eventnumber+1,200):
                pp.close()
                pp = PdfPages(directory + '_SavedEventsCommon_' + str(count) + '.pdf')
                count+=1
            pp.savefig(fig)
            print('{} out of {} saved!'.format(str(eventnumber), str(len(ind1))))
            fig.clear()
            plt.close(fig)
        pp.close()

def PlotRecursiveLPResults(RoughEventLocations, inp, directory, buffer, channel='i2'):
    pp = PdfPages(directory + '_' + channel + '_DetectedEventsFromLPFilter.pdf')
    a=1
    for i in RoughEventLocations['RoughEventLocations']:
        startp = np.uint64(i[0]-buffer*inp['samplerate'])
        endp = np.uint64(i[1]+buffer*inp['samplerate'])
        t = np.arange(startp, endp)
        t = t / inp['samplerate'] * 1e3
        fig = PlotEvent(t, inp[channel][startp:endp], channel=channel)
        pp.savefig(fig)
        print('{} out of {} saved!'.format(str(a), str(len(RoughEventLocations['RoughEventLocations']))))
        a+=1
        fig.clear()
        plt.close(fig)
    pp.close()

def SaveAllAxopatchEvents(AnalysisResults, directory, out, buffer, withFit = 1):
    # Plot All Common Events
    pp = PdfPages(directory + '_SavedEventsAxopatch.pdf')
    t = np.arange(0, len(out['i1']))
    t = t / out['samplerate'] * 1e3

    for eventnumber in range(AnalysisResults['i1']['NumberOfEvents']):
        parttoplot = np.arange(AnalysisResults['i1']['StartPoints'][eventnumber] - buffer,
                               AnalysisResults['i1']['EndPoints'][eventnumber] + buffer, 1, dtype=np.uint64)

        fit1 = np.concatenate([np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][eventnumber],
                               np.ones(AnalysisResults['i1']['EndPoints'][eventnumber] -
                                       AnalysisResults['i1']['StartPoints'][
                                           eventnumber]) * (
                                   AnalysisResults['i1']['LocalBaseline'][eventnumber] -
                                   AnalysisResults['i1']['DeltaI'][eventnumber]),
                               np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][eventnumber]])

        if withFit:
            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], fit1=fit1)
        else:
            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot])
        try:
            pp.savefig(fig)
        except:
            print('Problem at {} !'.format(str(eventnumber)))

        print('{} out of {} saved!'.format(str(eventnumber), str(AnalysisResults['i1']['NumberOfEvents'])))
        #print('Length i1: {}, Fit i1: {}'.format(len(out['i1'][parttoplot]), len(fit1)))
        fig.clear()
        plt.close(fig)
    pp.close()

def SetCusum2ARL0(deltax,sigmax,Arl0_2=1000,thresholdlevels=[1000,0.1,10]):
    h0min=thresholdlevels[0];
    h0max=thresholdlevels[1]; #detection threshold interval

    ARL0=2*ARL0_2; #for two-sided algo
    # Optimization on h "hopt"
    change=0;
    mu=deltax*(change-deltax/2)/sigmax**2;
    sigma=abs(deltax)/sigmax;
    mun=mu/sigma;
    def f(h):
        (exp(-2*mun*(h/sigma+1.166))-1+2*mun*(h/sigma+1.166))/(2*mun^2)-ARL0

    if(f(h0min)*f(h0max)<0):
        print('test')
        #hopt=fzero(f,[h0min h0max])
    elif(f(h0min) <0 and f(h0max) <0):
        hopt=h0max
    elif(f(h0min) >0 and f(h0max) >0):
        hopt=h0min

    hbook=sigmax*hopt/deltax;
    return hbook, hopt

def event_detection(RawSignal, CusumParameters, RoughEventLocations, cutoff):
    for i in range(len(RoughEventLocations)):
        CusumReferencedEndPoint = RoughEventLocations[i][1] + CusumParameters['BaselineLength'];
        CusumReferencedStartPoint = RoughEventLocations[i][0] - CusumParameters['BaselineLength'];

        if CusumReferencedEndPoint > len(RawSignal):
            CusumReferencedEndPoint = len(RawSignal)
        if CusumReferencedStartPoint < 0:
            CusumReferencedStartPoint = 0

        if RoughEventLocations[i][2] < CusumParameters['ImpulsionLimit'] * CusumParameters['SamplingFrequency']:
            trace = RawSignal[CusumReferencedStartPoint:CusumReferencedEndPoint]
             #[mc,krmv]=ImpulsionFitting(trace,BaselineLength, cutoff, SamplingFrequency);
             #EventType='Impulsion';
             #krmv=[BaselineLength, BaselineLength+RoughEventLocations(i,3)];
             #mc=[ones(1,BaselineLength+1)*mean(trace(1:BaselineLength)), ones(1,krmv(2)-krmv(1)-1)*min(trace), ones(1,BaselineLength-2)*mean(trace(krmv(2):end))];
        else:
            mc, kd, krmv = cusum(RawSignal[CusumReferencedStartPoint:CusumReferencedEndPoint],CusumParameters['delta'],CusumParameters['h'])
            EventType='Standard Event';

        Startpoint = CusumReferencedStartPoint + krmv[0]
        EndPoint = CusumReferencedStartPoint + krm[-1]

        BaselinePart = RawSignal[CusumReferencedStartPoint:CusumReferencedStartPoint]

def EredoxBefore14062018():
    Eredox = np.array([31.7, 82.9, 135, 185], dtype=float)
    Eredox = Eredox * 1e-3
    return (Eredox-Eredox[0])
