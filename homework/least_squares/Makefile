libraries = matrix linalg sfuns
library_path = ../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

Plot.svg: plot.gp out.txt
	gnuplot $<

out.txt: make_libraries main.bin
	./main.bin

main.bin: main.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f main.bin out.txt Plot.svg