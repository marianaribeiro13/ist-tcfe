# ------------------------------------------------------------
# type "make" command in the Linux terminal create report.pdf
#
# type "make clean" to delete all generated files
# ------------------------------------------------------------

all:
	make -C doc
	cp doc/report.pdf .

clean:
	make -C doc clean
	make -C sim clean
	make -C mat clean

.PHONY: all clean


