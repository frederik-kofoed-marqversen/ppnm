Out.txt: main.bin
	./main.bin > Out.txt

main.bin: main.rs
	rustc main.rs -o main.bin

.PHONEY: clean
clean:
	rm -f Out.txt  main.bin
