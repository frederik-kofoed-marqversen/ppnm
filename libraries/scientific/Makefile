OUT_DIR = target
lib_path = src/lib.rs

name = $(shell basename $(CURDIR))
library_files := $(shell find src -name '*.rs')
rlib_target = $(OUT_DIR)/lib$(name).rlib
test_target = $(OUT_DIR)/$(name).test

lib: $(rlib_target)

test:
	rustc $(lib_path) --test -o $(test_target)
	./$(test_target) monte_carlo
	rm $(test_target)

$(rlib_target): $(library_files)
	rustc $(lib_path) -O --crate-name $(name) --crate-type lib --out-dir $(OUT_DIR)

.PHONY: clean
clean:
	rm target/*