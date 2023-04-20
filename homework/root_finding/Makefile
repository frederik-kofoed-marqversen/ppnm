libraries = scientific matrix sfuns
library_path = ../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

run: make_libraries Out.txt Ground_state.svg Convergence.svg

Ground_state.svg: plot.gp Out.txt
	gnuplot $< > $@

Convergence.svg: plot_convergence.gp convergence.txt
	gnuplot $< > $@

convergence.txt: convergence.bin
	./$< > $@

Out.txt: main.bin
	./$< > $@

convergence.bin: convergence.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

main.bin: main.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f *.bin *.svg *.txt