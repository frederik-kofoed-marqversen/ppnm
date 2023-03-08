name = $(shell basename $(CURDIR))
lib_path = src/lib.rs
target_folder = ../target
library_files := $(wildcard src/*.rs)

lib: $(target_folder)/lib$(name).rlib

test: $(target_folder)/$(name).test
	./$<

$(target_folder)/lib$(name).rlib: $(library_files)
	rustc $(lib_path) --crate-name $(name) --crate-type lib --out-dir $(target_folder)

$(target_folder)/$(name).test: $(library_files)
	rustc $(lib_path) --test -o $@

.PHONEY: clean
clean:
	rm $(target_folder)/*$(name)*