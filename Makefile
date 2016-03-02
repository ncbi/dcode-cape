
LIB_DIR = src/libkm
KW_DIR = src/kweight
CD_DIR = src/cDescritors
FF_DIR = src/formatFasta

.PHONY: all
	
all: clean
	ln -s $(LIB_DIR)/includes/BedFactory.h includes/BedFactory.h
	ln -s $(LIB_DIR)/includes/FastaFactory.h includes/FastaFactory.h
	ln -s $(LIB_DIR)/includes/Global.h includes/Global.h
	ln -s $(LIB_DIR)/includes/KmersFactory.h includes/KmersFactory.h
	ln -s $(LIB_DIR)/includes/TimeUtils.h includes/TimeUtils.h
	ln -s $(LIB_DIR)/includes/SNPFactory.h includes/SNPFactory.h
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
	rm -rfv includes/*
	rm -rfv lib/*
	