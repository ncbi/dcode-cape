
LIB_DIR = src/libkm
KW_DIR = src/kweight
KM_DIR = src/kMerge
CD_DIR = src/cape
FF_DIR = src/formatFasta
SNP_DIR = src/snp2Fasta

.PHONY: all
	
all: clean
	$(MAKE) -C $(LIB_DIR)
	$(MAKE) -C $(KW_DIR)
	$(MAKE) -C $(KM_DIR)
	$(MAKE) -C $(CD_DIR)
	$(MAKE) -C $(FF_DIR)
	$(MAKE) -C $(SNP_DIR)
	echo "Done!!!"
	
clean:
	$(MAKE) -C $(SNP_DIR) clean
	$(MAKE) -C $(FF_DIR) clean
	$(MAKE) -C $(CD_DIR) clean
	$(MAKE) -C $(KM_DIR) clean
	$(MAKE) -C $(KW_DIR) clean
	$(MAKE) -C $(LIB_DIR) clean
	rm -rfv bin/*
	