libraries = scientific
library_path = ../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

run: make_libraries Out.txt Error_scaling.svg Stratified_samples.svg

Stratified_samples.svg: strat.gp strat.txt
	gnuplot $< > $@

Error_scaling.svg: plot.gp scaling.txt
	gnuplot $< > $@

strat.txt: strat.bin
	./$< > $@

scaling.txt: scaling.bin
	./$< > $@

Out.txt: main.bin
	./$< > $@

strat.bin: strat.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

scaling.bin: scaling.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

main.bin: main.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f *.bin *.svg *.txt