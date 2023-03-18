libraries = sfuns linalg matrix rand
library_path = ../../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

all: make_libraries Test.txt Gnulot_comparison.svg Plot.svg

Gnulot_comparison.svg: comparison.gp data.txt
	gnuplot $<

Plot.svg: plot.gp data.txt
	gnuplot $<

data.txt: main.bin
	./main.bin > $@

Test.txt: test.bin
	./test.bin > $@

main.bin: main.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

test.bin: test.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f *.bin *.svg *.txt