.cu.o:
	echo "NVCC_PATH: $NVCC_PATH"
	$(NVCC_PATH) -gencode=arch=compute_13,code=sm_13 -o $@ -c $<

.cu.lo:
	$(top_srcdir)/config/cudalt.py $@ $(NVCC_PATH) -gencode=arch=compute_13,code=sm_13 --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" -c $<

