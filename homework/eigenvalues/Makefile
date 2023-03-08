libraries = matrix rand sfuns
library_path = ../../libraries

Out.txt: main.bin
	./$< > $@

main.bin: main.rs
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
	rustc $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONEY: clean
clean:
	rm -f Out.txt  main.bin