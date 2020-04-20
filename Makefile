
SRCDIR=src
CPPFILE=reproject.cpp
ALGLIB = includes/alglibinternal.cpp includes/alglibmisc.cpp includes/ap.cpp includes/linalg.cpp includes/optimization.cpp includes/specialfunctions.cpp includes/statistics.cpp
ALGLIB2 = includes/alglibinternal.cpp includes/alglibmisc.cpp includes/ap.cpp includes/linalg.cpp

all: $(SRCDIR)/$(CPPFILE)
	g++ $(SRCDIR)/$(CPPFILE) -w -std=c++11 ${ALGLIB2}

clean:
	rm *~

