OUT_DIR=target
lib_path = src/lib.rs

name = $(shell basename $(CURDIR))
library_files := $(shell find src -name '*.rs')
rlib_target = $(OUT_DIR)/lib$(name).rlib

lib: $(rlib_target)

$(rlib_target): $(library_files) ../target/libnum_traits.rlib
	rustc --crate-name num_complex --edition=2018 $(lib_path) --crate-type lib --cfg 'feature="default"' --cfg 'feature="std"' --out-dir $(OUT_DIR) --extern num_traits=../target/libnum_traits.rlib

../target/libnum_traits.rlib:
	$(MAKE) -C ../ num_traits

.PHONY: clean
clean:
	rm target/*