#!/usr/bin/python

import sys,string
from subprocess import call
import shutil, os

file = open(sys.argv[1])
data = file.readlines()
file.close()

for line in data:
	aux = line.find('epsfile')
	if (aux == 0):
		filename = (line[9:-2])[:-4]
		break
print(filename)

call(['gnuplot' ,sys.argv[1]])

texfiledata=""" 
\\documentclass[prl,10pt]{revtex4-2}
\\usepackage[dvips]{graphicx}

\\begin{document}
\\thispagestyle{empty}

\\begin{figure}
\\input{"%(filename)s"}
\\end{figure} 

\\end{document}
"""%locals()

file = open("___aux.tex",'w')
file.write(texfiledata)
file.close()

call(['latex', '___aux.tex'])
call(['dvips', "___aux.dvi" ])
call(['ps2pdf', "___aux.ps" ])
call(['pdfcrop', "___aux.pdf" ])
shutil.move("___aux-crop.pdf", filename+".pdf")
os.system("rm ___*")


