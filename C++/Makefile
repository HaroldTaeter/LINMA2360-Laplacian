#
#  Makefile for Linux 
#
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#
#  A adapter en fonction des ordinateurs/environnements 
#  Compilateur, edition de liens, 
#
#
CC       = g++  
LD       = g++
CFLAGS   = -O3 -Dgraphic -g -Wall
LFLAGS   = -Wall -O3 -g
LIBS     = -lglfw -lm -lGL -lglut -lGLU -lpthread 
#
PROG     = myLaplacian.exe
LISTEOBJ = \
  laplacian.o main.o #glfem.o #graphic.o
# ATTENTION... aucun caractere apres le caractere de continuation "\"
#
# compilation
#
.c.o :
	$(CC) -c  $(CFLAGS) -o $@ $<
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# dependances
#
all        : $(PROG)
laplacian.o 	   : laplacian.cpp
main.o		   : main.cpp 
#glfem.o            : glfem.c delaunay.h glfem.h
#
# edition de lien
#
$(PROG) : $(LISTEOBJ)
	$(LD) -o $(PROG) $(LFLAGS) $(LISTEOBJ) $(LIBS)
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# effacement des fichiers intermediaires
#
clean :
	rm -vf $(PROG) $(LISTEOBJ) core a.out
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# ATTENTION... il faut une ligne vide a la fin du fichier.

