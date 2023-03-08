libraries = matrix sfuns rand complex

all: $(libraries)

.PHONY: $(libraries)
$(libraries): | target
	make -C $@ lib

target:
	mkdir target

.PHONY: clean
clean:
	rm target/*
