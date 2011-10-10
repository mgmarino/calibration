fit2.so: fit2.o DictOutput.o
	g++ -shared -o $@ `root-config --libs` -lMinuit $^

DictOutput.cpp: fit2.h LinkDef.h
	rootcint -f $@ -c $^

fit2.o: fit2.cpp
	g++ -fPIC -g `root-config --cflags` -c $^

DictOutput.o: DictOutput.cpp
	g++ -fPIC -g `root-config --cflags` -c $^

clean:
	rm -f fit2.so
	rm -f *.o
	rm -f DictOutput.*
