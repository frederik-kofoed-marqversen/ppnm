Out.txt: hello
	./hello > Out.txt

hello: hello.rs
	rustc hello.rs

.PHONEY: clean
clean:
	rm -f Out.txt  hello
