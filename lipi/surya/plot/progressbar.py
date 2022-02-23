import sys
import time

def ProgressBar(Total, Progress, BarLength=20, ProgressIcon="#", BarIcon="-"):
    '''
    Example: Bar=progressbar.ProgressBar(n,i);progressbar.ShowBar(Bar)
    '''
    try:
        # You can't have a progress bar with zero or negative length.
        if BarLength <1:
            BarLength = 20
        # Use status variable for going to the next line after progress completion.
        Status = ""
        # Calcuting progress between 0 and 1 for percentage.
        Progress = float(Progress) / float(Total)
        # Doing this conditions at final progressing.
        if Progress >= 1.:
            Progress = 1
            Status = "\r\n"    # Going to the next line
        # Calculating how many places should be filled
        Block = int(round(BarLength * Progress))
        # Show this
        Bar = "[{}] {:.0f}% {}".format(ProgressIcon * Block + BarIcon * (BarLength - Block), round(Progress * 100, 0), Status)
        return Bar
    except:
        return "ERROR"
        
def ShowBar(Bar):
    sys.stdout.write(Bar)
    sys.stdout.flush()
