# converts python (.py) to jupiter note book (.ipynb)
# howto: tojupiter.py program.py

from nbformat import v3, v4
import sys

content = open(sys.argv[1]).read()
outname = sys.argv[1].split('.')[0] + '.ipynb'
nb = v3.reads_py(content)
nb = v4.upgrade(nb)
nb_json = v4.writes(nb) + '\n'

outfile = open(outname, 'w')
outfile.write(nb_json)
outfile.close()
