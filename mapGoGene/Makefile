SRC =  AnnotMgr.C AnnotTerm.C Framework.C GoTerm.C GoTermManager.C
LIBPATH =/mnt/ws/sysbio/roygroup/shared/thirdparty/gsl_install/lib
INCLPATH1 = /mnt/ws/sysbio/roygroup/shared/thirdparty/gsl_install/include
LFLAG = -lgsl -lgslcblas 

CC=g++
CFLAGS = -g

mapGoGene: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o mapGoGene
	cp mapGoGene ../../
clean: 
	rm mapGoGene
