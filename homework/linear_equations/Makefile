libraries = matrix rand sfuns

library_path = ../../libraries
rlib_path := $(library_path)/target
find_rlibs = $(rlib_path)/lib$(lib).rlib

Out.txt: main.bin
	./main.bin > $@

main.bin: main.rs $(foreach lib,$(libraries),$(find_rlibs))
	rustc main.rs $(foreach lib,$(libraries),--extern $(lib)=$(find_rlibs)) -o $@

$(rlib_path)/lib%.rlib :: $(library_path)/%.rs
	rustc $< --crate-type=lib -o $@

.PHONEY: clean
clean:
	rm -f Out.txt  main.bin
