AUTOMAKE_OPTIONS = foreign
SUBDIRS = lib bin test

scan:
	script/scan.py lib > lib/Makefile.am
	script/scan.py bin > bin/Makefile.am
	echo > test/Makefile.am

scan_test:
	script/scan.py lib > lib/Makefile.am
	script/scan.py bin > bin/Makefile.am
	script/scan.py test > test/Makefile.am
