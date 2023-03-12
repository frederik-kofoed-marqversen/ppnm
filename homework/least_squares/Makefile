libraries = matrix linalg sfuns
library_path = ../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

plot: Out.txt plot.gp
	gnuplot plot.gp

Out.txt: make_libraries main.bin
	./main.bin

main.bin: main.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f main.bin Out.txt plot.svg