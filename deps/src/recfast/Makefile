LIBRECFAST = ../../usr/lib/librecfast.so

default:
	gfortran -ffixed-line-length-none -fPIC -shared -g -O0 recfast.for -o ${LIBRECFAST}
	
clean:
	rm -f ${LIBRECFAST}
