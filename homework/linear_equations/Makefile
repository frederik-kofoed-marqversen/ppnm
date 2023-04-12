libraries = matrix sfuns scientific
library_path = ../../libraries

all: Out/test.txt Out/plot.svg

Out/plot.svg: plot.gp Out/times.data
	gnuplot $< > $@

Out/test.txt: main.bin
	./$< -test > $@

Out/times.data : main.bin
	>$@
	for N in $$(seq 200 200 2000); do \
		echo "Solving with size $$N <= 2000"; \
		time --format "$$N\t%e" --output $@ --append \
		./$< -size:$$N;\
	done

main.bin: main.rs
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
	rustc $< -O $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONEY: clean
clean:
	rm -f main.bin Out/*