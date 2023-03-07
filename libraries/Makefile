libraries = matrix rand sfuns
target_pattern = target/lib%.rlib
targets := $(libraries:%=$(target_pattern))

all: $(libraries)

$(libraries) : % : $(target_pattern)

$(target_pattern) :: %.rs
	rustc $< --crate-type=lib -o $@

.PHONEY: clean
clean:
	rm target/*
