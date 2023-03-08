OUT_DIR=target
lib_path = src/lib.rs

name = $(shell basename $(CURDIR))
library_files := $(shell find src -name '*.rs')
rlib_target = $(OUT_DIR)/lib$(name).rlib

lib: $(rlib_target)

$(rlib_target): $(library_files)
	rustc --crate-name num_traits $(lib_path) --crate-type lib --cfg 'feature="std"' --out-dir $(OUT_DIR)

.PHONY: clean
clean:
	rm target/*