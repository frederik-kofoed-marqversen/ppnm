libraries = sfuns
library_path = ../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

all: make_libraries Erf.svg Out.txt Tests.txt

Tests.txt: test.bin
	./$< > $@
	
test.bin: main.rs
	rustc $< --test -o $@ $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib)

Erf.svg: erf.gp erf.data
	gnuplot $< > $@

erf.data Out.txt: main.bin
	./$<

main.bin: main.rs
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f *.bin *.svg *.data *.txt