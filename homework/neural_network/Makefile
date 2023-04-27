libraries = matrix sfuns scientific
library_path = ../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

run: make_libraries Curve_fit.svg Diff_eq.svg

Curve_fit.svg Diff_eq.svg: plot.gp data.txt
	gnuplot $<

data.txt: main.bin
	./$< > $@

main.bin: main.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f *.bin *.svg *.txt