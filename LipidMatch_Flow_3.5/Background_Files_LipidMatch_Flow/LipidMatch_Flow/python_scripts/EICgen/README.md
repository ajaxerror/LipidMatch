# EICgen
This script to generate EIC CSV files must be compiled to run 
with the LipidMatch project

### creating EICgen executable

ON WINDOWS: run

`pip install -r requirements.txt`
`pyinstaller EIC_gen.py`.

Preferably, you would do this inside a virtual environment and test the script before compiling.
The pyinstaller packager will build the dists directory and EICgen bundle within it. The bundle
will contain an executable that can be called by LipidMatch.r.
