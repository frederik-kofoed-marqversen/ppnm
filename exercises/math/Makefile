Out.txt: main.bin
	./main.bin > $@

main.bin: main.rs libsfuns.rlib
	rustc main.rs --extern sfuns=libsfuns.rlib -o $@

libsfuns.rlib: sfuns.rs
	rustc --crate-type=lib $< -o $@

.PHONEY: clean
clean:
	rm -f Out.txt  main.bin libsfuns.rlib
