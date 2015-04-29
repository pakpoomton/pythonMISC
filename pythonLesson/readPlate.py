import csv
import numpy as np
import re
import matplotlib.pyplot as plt

### SETUP INPUT PARAMETERS ACCORDING TO EXPERIMENT 
# Name of cvs file to read from
fileName = '042815plate.csv'
# Well on the 96 well plate we use. Here we use well number 1,2,...,7
selectedWell = np.array(range(1,8))
## Specify plot display. Note that we read from 7 well for this experiment but
## The first one is blank (LB with cell) which we use for background subtraction
## and won't be plotted. Thus, the first element of colorList, lst and leb won't
## be used for anything
# list of plot color corresponding to data from each well in (red, green, blue)
colorList = [(0,0.5,0), (0,0,1), (0,0,1), (1,0,0), (1,0,1), (1,0,1), (0,0.7,0)]
# line style of plot corresponding to data from each well
lst = ['-','--', '-', '-', '--','-','-']
# label on each plot to be put in legend 
leb = ['B','117','119','179','117/179','119/179','SAR20']
## set the upper bound of Y-axis for each plot. For now, we will manually adjust it
## according to dataset
yMax = [80000, 250000, 5000, 3.5]
## set Y-axis label according to channelName
channelName = ['YFP', 'CFP', 'RFP', 'OD600']


###----------------------------------------------------
## ACQUIRE AND ORGANISE DATA
# read from cvs to array of list
array = list(csv.reader(open(fileName, 'r')))
# find out which row in CSV file the data table start. 
sampleName = [row[2] for row in array] # get the third column
# The table start at the row below the row whose third column has a word 'Content')
sampleStartIndex = sampleName.index('Content')+1 





def getChannelInd(txt):
    ## Extract measurement channel index from data table header
    # data header format has measurement channel index inside '()'. 
    # For example, one of the header of data column is
    # 'Raw Data (600, UVlens 4)\n58 - 11 h 4 min '
    # number '4' inside the parenthesis indicates that this column belongs
    # to 4th measurement channel.
    # To group raw data by channel, we first define regular expression
    # note: need \ as ) is a special character
    closeParenthesis = re.compile('\)')
    m = closeParenthesis.search(txt)
    # find where ')' locate
    closeParenthesisRange = m.span()
    # get channel index number    
    channelInd = int(txt[closeParenthesisRange[0] - 1])
    return channelInd

def getTime(txt):
    ## Extract measurement time from data table header
    # data header format has time behind '-' in hour and min
    #  (see ex in getChannel comment)
    # to extract time in hour and minute, we first define RE pattern
    hourPattern = re.compile('-.*h')
    mHr = hourPattern.search(txt)
    locHr = mHr.span()
    hours = float(txt[locHr[0]+1:locHr[1]-1])
    minPattern = re.compile('h.*m')
    mMin = minPattern.search(txt)
    if mMin != None:
        locMin = mMin.span()
        mins = float(txt[locMin[0]+1:locMin[1]-1])
    else:
        mins = 0.
    totalMins = hours*60+mins
    return totalMins

## Get the array of measurement time
# extract the header text of data table 
headerRow = array[sampleStartIndex-1]
dataHeaderRow = headerRow[3:]

# Apply getTime to every entry in header row
measuredTime = map(getTime, dataHeaderRow)
tMax = max(measuredTime) # last time point to be used in the plot

## Get the starting location of each measurement channel 
## Find which rows have data from the well we want
selectedIndex = sampleStartIndex + selectedWell

## Find which columns has the first & last time point of ech channel measure 
# apply getChennelInd to every entry in header row
channelIndList = map(getChannelInd, dataHeaderRow)
# find the beginning column of each channel
channelDiff =  np.diff(channelIndList, n=1, axis=0)
channelChange =  channelDiff.ravel().nonzero()
channelChange = np.array([0]+list(channelChange[0]+1)+[len(dataHeaderRow)])

## ----------------------------------
## GENERATING PLOT
# get measurement results from blank well (LB without cell)
dataBlank = array[selectedIndex[0]]
measuredBlank = map(float, dataBlank[3:])

# Iterate through measurement from different wells and plot a line for each well
for n in range(1,len(selectedIndex)):
    
    well_n = selectedIndex[n]
    dataRow = array[well_n] # measurement data from n_th well 
    measuredData = np.array(map(float, dataRow[3:]))-np.array(measuredBlank)# subtract blank
    channel1Data = measuredData[channelStartInd[0]:channelStartInd[1]]
    channel2Data = measuredData[channelStartInd[1]:channelStartInd[2]]
    channel3Data = measuredData[channelStartInd[2]:channelStartInd[3]]
    channel4Data = measuredData[channelStartInd[3]:channelStartInd[4]]
    channel1Time = measuredTime[channelStartInd[0]:channelStartInd[1]]
    channel2Time = measuredTime[channelStartInd[1]:channelStartInd[2]]
    channel3Time = measuredTime[channelStartInd[2]:channelStartInd[3]]
    channel4Time = measuredTime[channelStartInd[3]:channelStartInd[4]]
    
    ax = plt.subplot(411)
    plt.plot(channel1Time, channel1Data , color = colorList[n], label=leb[n], linewidth=2.0, linestyle=lst[n])
    plt.ylabel(channelName[0])
    plt.xticks([])
    plt.axis([0, tMax, -5000, yMax[0]])
    
    plt.subplot(412)
    plt.plot(channel2Time, channel2Data , color = colorList[n], linewidth=2.0, linestyle=lst[n])
    plt.ylabel(channelName[1])
    plt.xticks([])
    plt.axis([0, tMax, -10000, yMax[1]])
      
    plt.subplot(413)
    plt.plot(channel3Time, channel3Data , color = colorList[n], linewidth=2.0, linestyle=lst[n])
    plt.ylabel(channelName[2])
    plt.xticks([])
    plt.axis([0, tMax, -200, yMax[2]])
      
    plt.subplot(414)
    plt.plot(channel4Time, channel4Data , color = colorList[n], linewidth=2.0, linestyle=lst[n])
    plt.ylabel(channelName[3])
    plt.xlabel('time (min)')
    plt.axis([0, tMax, -0.3, yMax[3]])


ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.5),ncol=3)
plt.show()


