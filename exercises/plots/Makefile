libraries = sfuns

library_path = ../../libraries
rlib_path := $(library_path)/target
find_rlibs = $(rlib_path)/lib$(lib).rlib

plots: target/gamma.svg target/lngamma.svg

# Do the plots and save to /target/ folder
target/gamma.svg target/lngamma.svg: plot.gp target/gamma.data target/factorials.data target/lngamma.data
	gnuplot $<
	echo .svg files has been put in the /target/ folder

# Run binary to optain the .data files
target/gamma.data target/lngamma.data: target/main.bin
	./$<

# Compile the rust binary
target/main.bin: main.rs $(foreach lib,$(libraries),$(find_rlibs))
	rustc main.rs $(foreach lib,$(libraries),--extern $(lib)=$(find_rlibs)) -o $@

$(rlib_path)/lib%.rlib :: $(library_path)/%.rs
	rustc $< --crate-type=lib -o $@

.PHONEY: clean
clean:
	rm target/*