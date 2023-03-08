libraries = matrix sfuns rand complex matrix_v0

all: $(libraries)

.PHONY: $(libraries)
$(libraries): | target
	$(MAKE) -C $@ lib OUT_DIR=../target > /dev/null

target:
	mkdir target

.PHONY: clean
clean:
	rm target/*
