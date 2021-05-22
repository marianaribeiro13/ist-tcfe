all: nspice.log pspice.log
	$(eval NGSFIGS=$(shell grep _FIG nspice.log| sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval NGSFIG=$(shell grep _FIG pspice.log| sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGP=$(addsuffix .ps, $(NGSFIG)))
	$(foreach i, $(NGSFIGP), ps2pdf $i;)

nspice.log:../sim/t4.net
	ngspice -b $< -o $@

pspice.log: ../sim/t4aux.net
	ngspice -b $< -o $@

clean:
	@rm -f nspice.log pspice.log *.tex *.pdf *.ps *~

.PHONY: all clean
