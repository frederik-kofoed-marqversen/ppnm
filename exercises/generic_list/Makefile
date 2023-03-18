run: main.bin data.txt
	./$< -i:data.txt

main.bin: main.rs
	rustc $< -o $@

.PHONY: clean
clean:
	rm -f main.bin