libraries = sfuns

library_path = ../../libraries
rlib_path := $(library_path)/target
find_rlibs = $(rlib_path)/lib$(lib).rlib

plots: out/gamma.svg out/lngamma.svg
	$(info .svg files has been put in the out/ folder)

# Do the plots and save to /out/ folder
out/gamma.svg out/lngamma.svg: plot.gp out/gamma.data out/factorials.data out/lngamma.data
	gnuplot $<

# Run binary to optain the .data files
out/gamma.data out/lngamma.data: main.bin
	./$<

# Compile the rust binary
main.bin: main.rs $(foreach lib,$(libraries),$(find_rlibs))
	rustc main.rs $(foreach lib,$(libraries),--extern $(lib)=$(find_rlibs)) -o $@

$(rlib_path)/lib%.rlib :: $(library_path)/%.rs
	rustc $< --crate-type=lib -o $@

.PHONEY: clean
clean:
	rm out/* main.bin