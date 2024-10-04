PYTHON = python3
FLAKE8 = python3 -m flake8

# for now, just cupcake imports
cupcake_src = sqanti3/cupcake/*/*.py
cupcake_ignore = C901,E101,E111,E126,E127,E122,E125,E128,E202,E203,E221,E225,E226,E231,E261,E275,E262,E265,E301,E302,E303,E305,E502,E701,E704,E722,W191,W504

all: lint

lint:
	${FLAKE8} --color=never --ignore=${cupcake_ignore} ${cupcake_src}

