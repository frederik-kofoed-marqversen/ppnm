test: main.bin input.txt
	./main.bin -numbers:1,2,3,4,5 -i:input.txt -o:output.txt < In.txt > Out.txt

main.bin: main.rs
	rustc main.rs -o main.bin

.PHONEY: clean
clean:
	rm -f Out.txt main.bin output.txt
