libraries = matrix rand sfuns
library_path = ../../libraries
rlibs := $(libraries:%=$(library_path)/target/lib%.rlib)

drs = 0.03 0.04 0.06 0.08 0.1 0.15 0.2 0.3 0.4
r_maxs = 10 9 8 7 6 5 4 3 2
sizes = 350 325 300 275 250 225 200 175 150 125 100 75 50 25

all: tests hydrogen_plots

# TEST METHOD
tests: Out/test.txt Out/times.svg

Out/times.svg: plot_tests.gp times.data
	gnuplot $< > $@

Out/test.txt: test.bin
	./test.bin -print:true -size:5 > $@

times.data: test.bin
	> $@
	$(MAKE) $(sizes:%=size-%)
	echo "\n" >> $@
	$(MAKE) $(sizes:%=size2-%)

size-%: test.bin
	time --format "$(subst size-,,$@) %e" --output times.data --append ./$< -size:$(subst size-,,$@)

size2-%: test.bin
	time --format "$(subst size2-,,$@) %e" --output times.data --append ./$< -size:$(subst size2-,,$@) -opt:true

# MAKE PLOTS OF HYDROGEN
hydrogen_plots: Out/dr_conv.svg Out/r_max_conv.svg Out/states.svg

Out/dr_conv.svg Out/r_max_conv.svg Out/states.svg: plot.gp dr.data r_max.data states.data
	gnuplot $<

dr.data: main.bin
	> $@
	$(MAKE) $(drs:%=dr-%)

r_max.data: main.bin
	> $@
	$(MAKE) $(r_maxs:%=r_max-%)

states.data: main.bin
	./$< -r_max:10 -dr:0.05 -N:3 > $@

dr-%: main.bin
	./$< -r_max:10 -dr:$(subst dr-,,$@) >> dr.data

r_max-%: main.bin
	./$< -r_max:$(subst r_max-,,$@) -dr:0.05 >> r_max.data

# MAKE .bins
main.bin: main.rs | make_libraries
	rustc $< -O $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

test.bin: test.rs | make_libraries
	rustc $< -O $(foreach lib,$(libraries),--extern $(lib)=$(library_path)/target/lib$(lib).rlib) -o $@

# CLEAN AND .rlibs
.PHONEY: clean make_libraries
make_libraries:
	$(MAKE) -C $(library_path) $(libraries) > /dev/null
clean:
	rm -f *.bin *.data Out/*