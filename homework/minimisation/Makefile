libraries = matrix
library_path = ../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

run: make_libraries Out.txt Higgs.svg

Higgs.svg: plot.gp higgs.data Out.txt
	gnuplot $< > $@

Out.txt: main.bin
	./$< < higgs.data > $@

main.bin: main.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f *.bin *.svg *.txt