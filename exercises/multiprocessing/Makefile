N = 100000000
TIME = time --portability --append --output $@

Out.txt : main.bin Makefile
	> $@
	$(TIME) ./$< -terms:$N -threads:1 >> $@
	$(TIME) ./$< -terms:$N -threads:2 >> $@
	$(TIME) ./$< -terms:$N -threads:3 >> $@
	$(TIME) ./$< -terms:$N -threads:4 >> $@

main.bin: main.rs
	rustc $< -o $@

.PHONEY: clean
clean:
	rm -f Out.txt  main.bin
