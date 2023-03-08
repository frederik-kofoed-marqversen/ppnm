libraries = matrix sfuns rand num_complex num_traits

all: $(libraries)

.PHONY: $(libraries)
$(libraries): | target
	$(MAKE) -C $@ lib OUT_DIR=../target > /dev/null

target:
	mkdir target

.PHONY: clean
clean:
	rm target/*
