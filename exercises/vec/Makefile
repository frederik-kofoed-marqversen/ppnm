Out.txt: main.bin
	./main.bin > $@

main.bin: main.rs libvec3d.rlib
	rustc main.rs --extern vec3d=libvec3d.rlib -o $@

libvec3d.rlib: vec3d.rs
	rustc --crate-type=lib $< -o $@

.PHONEY: clean
clean:
	rm -f Out.txt  main.bin libvec3d.rlib
