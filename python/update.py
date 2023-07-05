import sys
#Method borrowed from Cosma Fulga who in turn obtained from the link referred below.

def update_progress(progress, decimalpoints=0):
    """ Make an interactive progress bar as described on:
    https://stackoverflow.com/questions/3160699/python-progress-bar
    """
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + 
            "-"*(barLength-block), round(progress*100, decimalpoints), status)
    sys.stdout.write(text)
    sys.stdout.flush()


