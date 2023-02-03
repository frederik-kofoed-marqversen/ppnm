Out.txt: main
	./main > Out.txt

main: main.rs
	rustc main.rs

.PHONEY: clean test-shell
clean:
	rm -f Out.txt  main
