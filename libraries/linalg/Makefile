OUT_DIR = target
lib_path = src/lib.rs
dependencies = matrix

name = $(shell basename $(CURDIR))
library_files := $(shell find src -name '*.rs') $(dependencies:%=../target/lib%.rlib)
rlib_target = $(OUT_DIR)/lib$(name).rlib
test_target = $(OUT_DIR)/$(name).test

lib: make_dependencies $(rlib_target)

test: make_dependencies
	rustc $(lib_path) --test -o $(test_target)
	./$(test_target)
	rm $(test_target)

$(rlib_target): $(library_files)
	rustc $(lib_path) --crate-name $(name) --crate-type lib --out-dir $(OUT_DIR) \
	$(foreach lib,$(dependencies),--extern $(lib)=../target/lib$(lib).rlib)

.PHONY: clean make_dependencies
make_dependencies:
	make -C ../ $(dependencies)
clean:
	rm target/*