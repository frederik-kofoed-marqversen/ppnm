libraries = matrix sfuns num_complex num_traits scientific

all: $(libraries)

.PHONY: $(libraries)
$(libraries): | target
	$(MAKE) -C $@ lib OUT_DIR=../target

target:
	mkdir target

.PHONY: clean
clean:
	rm target/*
