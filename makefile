
#------------------------------------------------------------------------------

SOURCE=yeeScheme.cc
MYPROGRAM=yee
#MYINCLUDES=/home/scale/g++Projects/gLib/
MYINCLUDES=/usr/include/boost/
#MYINCLUDES=/usr/include/boost/

#MYLIBRARIES=fltk
CC=g++

#------------------------------------------------------------------------------

all: $(MYPROGRAM)


$(MYPROGRAM): $(SOURCE)

	$(CC) -I$(MYINCLUDES) $(SOURCE) -o$(MYPROGRAM) 

clean:

	rm -f $(MYPROGRAM)

