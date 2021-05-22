all: ngspice.log 
	$(eval NGSFIGS=$(shell grep _FIG nspice.log| sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' spice.log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tabs.tex;)
	$(eval NGSFIG=$(shell grep _FIG pspice.log| sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGP=$(addsuffix .ps, $(NGSFIG)))
	$(foreach i, $(NGSFIGP), ps2pdf $i;)

ngspice.log: ../sim/t4.net 
	ngspice -b $< -o $@ 


clean:
	@rm -f ngspice.log *.tex *.pdf *.ps *~

.PHONY: all clean
