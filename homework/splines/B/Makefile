libraries = sfuns
library_path = ../../../libraries
rlib_files = $(libraries:%=$(library_path)/target/lib%.rlib)

Out.txt: main.bin
	./$< > $@

main.bin: main.rs | make_libraries
	rustc -A dead_code $< $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

.PHONY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f main.bin Out.txt