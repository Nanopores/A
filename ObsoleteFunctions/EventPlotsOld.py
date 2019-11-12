from bokeh.models import Legend,LegendItem,Range1d
from bokeh.plotting import figure, output_file, show, output_notebook
from bokeh.layouts import row, column, gridplot


def PlotGTauVoltageOld(eventClass, xLim=None, yLim=None, showCurrent=False, voltageLimits = None):
    bokeh.io.output_notebook(INLINE)

    #sort the voltages
    # categorize events in three types
    # CUSUM fitted events
    voltagesList = eventClass.GetAllVoltages()
    if voltageLimits is not None:
        voltagesList = [x for x in voltagesList if x >= voltageLimits[0] and x <= voltageLimits[1]]

    #define of variables
    TOOLS = "box_zoom,pan,wheel_zoom,reset"
    colors = sns.color_palette(n_colors=len(voltagesList))
    linecolors = sns.color_palette("muted", n_colors=len(voltagesList))

    assert(len(voltagesList)<=len(colors))
    backgroundColor = "#fafafa"
    bins = 50

    output_notebook()
    p = figure(plot_width=500, plot_height=500, min_border=10, min_border_left=50, tools=TOOLS,
               x_axis_location=None, y_axis_location=None,title="Linked Histograms")

    p.background_fill_color = "#fafafa"

    #select min max
    events = [event for event in eventClass.events if showCurrent or event.voltage is not 0]

    allTau, allYVals = extractytau(events, showCurrent)

    #Set limits
    if xLim is None:
        taurange = np.linspace(min(allTau), max(allTau), num=bins)
    else:
        assert(len(xLim) is 2)
        p.x_range=Range1d(xLim[0], xLim[1])
        taurange = np.linspace(xLim[0], xLim[1], num=bins)

    if yLim is None:
        yValsrange = np.linspace(min(allYVals), max(allYVals), num=bins)
    else:
        assert(len(yLim) is 2)
        p.y_range=Range1d(yLim[0], yLim[1])
        yValsrange = np.linspace(yLim[0], yLim[1], num=bins)

    #initialize values
    zerostau = np.zeros(len(taurange) - 1)
    previoustau = np.zeros(len(taurange) - 1)
    zerosy = np.zeros(len(yValsrange)-1)
    previousy = np.zeros(len(yValsrange) - 1)


    #Define histogram bottom
    ph = figure(plot_width=p.plot_width, plot_height=100, x_range=p.x_range,
                y_range=(0, 1), min_border=10, min_border_left=50, y_axis_location="right")

    ph.background_fill_color = backgroundColor
    ph.xgrid.grid_line_color = None
    ph.yaxis.major_label_orientation = np.pi / 4

    #Define the second histogram
    pv = figure(plot_width=300, plot_height=p.plot_height, x_range=(0, 1),
                y_range=p.y_range, min_border=10, y_axis_location="right")

    pv.ygrid.grid_line_color = None
    pv.xaxis.major_label_orientation = np.pi / 4
    pv.background_fill_color = backgroundColor

    legend_it = []

    i = 0
    #for i in range(len(voltagesList)):
    for voltage in voltagesList:
        color = "#%02x%02x%02x" % (int(colors[i][0]*255), int(colors[i][1]*255), int(colors[i][2]*255))
        linecolor = "#%02x%02x%02x" % (int(linecolors[i][0]*255), int(linecolors[i][1]*255), int(linecolors[i][2]*255))
        i+=1
        #for voltage, color, linecolor in zip(voltagesList, colors, linecolors):
        selectEvents = eventClass.GetEventsforVoltages(voltage)
        tau, yVals = extractytau(selectEvents, showCurrent)

        c = p.scatter(tau, yVals, size=3, line_width=2, color=color, alpha=0.8)

        hhist, hedges = np.histogram(tau, bins=taurange)

        hh = ph.quad(bottom=zerostau +previoustau, left=hedges[:-1], right=hedges[1:], top=hhist+previoustau,
                fill_color = color,line_color =  linecolor)

        previoustau +=hhist

        vhist, vedges = np.histogram(yVals, bins=yValsrange)
        vh = pv.quad(left=zerosy + previousy, bottom=vedges[:-1], top=vedges[1:], right=vhist+previousy,
                color=color, line_color= linecolor)
        previousy +=vhist

        #Sharing legend is broke
        #legend_it.append((name, [c, hh, vh]))
        name = str(Volt.format_data(voltage))
        legend_it.append(LegendItem(label=name, renderers=[vh]))

    #p.legend.location =  "top_right"
    #legend = Legend(items=legend_it, location=(0, -30))
    legend = Legend(items=legend_it)
    pv.add_layout(legend, 'right')

    #pv.legend.click_policy = "hide"

    ph.y_range.end = max(previoustau)*1.1
    pv.x_range.end = max(previousy) * 1.1

    fig = gridplot([[p, pv],[ph, None]], toolbar_location="left")

    # show the results
    show(fig)
