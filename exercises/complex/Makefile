libraries = num_complex num_traits
library_path = ../../libraries

Out.txt: main.bin
	./$< > $@

main.bin: main.rs
	$(MAKE) -C $(library_path) $(libraries)
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONEY: clean
clean:
	rm -f Out.txt  main.bin