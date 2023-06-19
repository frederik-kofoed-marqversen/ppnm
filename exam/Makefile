libraries = scientific matrix sfuns
library_path = ../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

all: make_libraries Test.txt Plot.svg

Plot.svg: plot.gp times.txt
	gnuplot $< > $@

Test.txt: main.bin
	./$< -test > $@

times.txt: main.bin
	>$@
	for N in $$(seq 100 50 1000); do \
		time --format "$$N\t%e" --output $@ --append \
		./$< -size:$$N -method:hessenberg &\
	done;\
	echo "\nTiming Hessenberg decomposition computation";\
	wait
	echo "" >> $@
	echo "" >> $@
	for N in $$(seq 100 50 1000); do \
		time --format "$$N\t%e" --output $@ --append \
		./$< -size:$$N -method:qr &\
	done;\
	echo "\nTiming QR factorisation computation";\
	wait

main.bin: main.rs $(rlib_files)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -C opt-level=0 -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f *.bin *.svg *.txt