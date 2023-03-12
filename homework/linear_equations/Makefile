libraries = matrix rand sfuns
library_path = ../../libraries

all: out/test.txt out/times.data out/plot.svg

out/plot.svg: plot.gp
	gnuplot plot.gp > out/plot.svg

plot.gp: out/times.data make_gp_file.sh
	./make_gp_file.sh $$(tail -n 1 $< | grep -oE '[^ ]+$$')

out/test.txt: main.bin
	./$< -test > $@

out/times.data : main.bin
	>$@
	for N in $$(seq 200 200 2000); do \
		echo "Solving with size $$N <= 2000"; \
		time --format "$$N %e" --output $@ --append \
		./$< -size:$$N;\
	done

main.bin: main.rs
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
	rustc $< -O $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONEY: clean
clean:
	rm -f main.bin out/*