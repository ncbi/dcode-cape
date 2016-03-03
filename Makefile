
LIB_DIR = src/libkm
KW_DIR = src/kweight
CD_DIR = src/cape
FF_DIR = src/formatFasta

.PHONY: all
	
all: clean
	$(MAKE) -C $(LIB_DIR)
	$(MAKE) -C $(KW_DIR)
	$(MAKE) -C $(CD_DIR)
	$(MAKE) -C $(FF_DIR)
	echo "Done!!!"
	
clean:
	$(MAKE) -C $(FF_DIR) clean
	$(MAKE) -C $(CD_DIR) clean
	$(MAKE) -C $(KW_DIR) clean
	$(MAKE) -C $(LIB_DIR) clean
	rm -rfv bin/*
	