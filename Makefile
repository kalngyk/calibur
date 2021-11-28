COMPILER = g++
INCL_DIR =  
HEADERS = InitCluster.h rmsd.h SimpPDB.h jacobi.h cubic.h
LIBRARY = 
SOURCES =
#CFLAGS= -O2 -D_USE_FAST_RMSD_ -D_SHOW_PERCENTAGE_COMPLETE_ -D_LARGE_DECOY_SET_
CFLAGS= -O2 -D_USE_FAST_RMSD_ -D_LARGE_DECOY_SET_
CONCERTLIBDIR = 

OBJECTS =  jacobi.o cubic.o rmsd.o SimpPDB.o PreloadedPDB.o main.o
SRC_PACKAGE_FILES = *.h *.cc Makefile README HISTORY

#a: PreloadedPDB.o SimpPDB.o
#	$(COMPILER) PreloadedPDB.o SimpPDB.o

all: calibur calibur-lite

calibur: $(OBJECTS) $(SOURCES) $(HEADERS) $(INCL_DIR) InitCluster.cc
	$(COMPILER) $(CFLAGS) $(INCL_DIR) -c  InitCluster.cc
	$(COMPILER) -o  $@ $(CFLAGS) $(OBJECTS) InitCluster.o $(LIBRARY) -lm 

calibur-lite: $(OBJECTS) $(SOURCES) $(HEADERS) $(INCL_DIR) InitCluster.cc
	$(COMPILER) $(CFLAGS) $(INCL_DIR) -D_ADD_LITE_MODE_ -c  InitCluster.cc
	$(COMPILER) -o  $@ $(CFLAGS) $(OBJECTS) InitCluster.o $(LIBRARY) -lm 

package:
	cp windows/calibur/bin/Debug/calibur.exe .
	zip -r calibur.win32.zip calibur.exe readme.txt
	rm calibur.exe
	rm -r calibur; :
	mkdir calibur
	cp $(SRC_PACKAGE_FILES) calibur
	tar zcf calibur.tar.gz calibur
	rm -r calibur
	cp $(SRC_PACKAGE_FILES) windows/calibur
	zip calibur.zip `cat windows/files`

%.o:%.cc
	$(COMPILER) $(CFLAGS) $(INCL_DIR) -c  $<

clean:
	rm *.o *.exe

