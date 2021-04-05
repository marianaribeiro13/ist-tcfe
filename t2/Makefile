# ------------------------------------------------------------
# type "make" command in the Linux terminal create report.pdf
#
# type "make clean" to delete all generated files
# ------------------------------------------------------------

all:
	@python t2_datagen.py > data.txt
	make -C mat
	make -C sim
	make -C doc
	cp doc/report.pdf . 

clean:
	@rm -f data.txt
	make -C doc clean
	make -C sim clean
	make -C mat clean
.PHONY: all clean


