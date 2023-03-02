Out.txt: target/main.bin
	./$< > $@

target/main.bin: main.rs target/libnum_traits.rlib target/libnum_complex.rlib
	rustc $< --edition=2021 --extern num_traits=target/libnum_traits.rlib --extern num_complex=target/libnum_complex.rlib -o $@

target/libnum_complex.rlib: ../../libraries/complex/num-complex target/libnum_traits.rlib
	rustc --crate-name num_complex --edition=2018 $</lib.rs --crate-type lib --cfg 'feature="default"' --cfg 'feature="std"' --out-dir target --extern num_traits=target/libnum_traits.rlib

target/libnum_traits.rlib: ../../libraries/complex/num-traits
	rustc --crate-name num_traits $</lib.rs --crate-type lib --cfg 'feature="std"' --out-dir target

.PHONEY: clean
clean:
	rm -f Out.txt  main.bin
	rm target/*
