libraries = matrix sfuns rand num_complex num_traits linalg

all: $(libraries)

.PHONY: $(libraries)
$(libraries): | target
	$(MAKE) -C $@ lib OUT_DIR=../target

target:
	mkdir target

.PHONY: clean
clean:
	rm target/*
