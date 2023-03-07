library_files := $(wildcard src/*/*.rs)
lib_path = src/lib.rs

test: target/matrix.test
	./$<

lib: target/libmatrix.rlib

target/matrix.test: $(library_files)
	rustc $(lib_path) --test -o $@

target/libmatrix.rlib: $(library_files)
	rustc $(lib_path) --crate-type lib -o $@

.PHONEY: clean
clean:
	rm target/*