# uses tojupyter.py to convert a .py to a .ipynb file
# howto: bash open_jupyter.sh filename.py

python tojupyter.py $1
name=${1:0:${#1}-3}.ipynb
jupyter notebook $name
