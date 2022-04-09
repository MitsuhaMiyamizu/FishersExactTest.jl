CXX = clang++
CXXFLAGS = -std=c++17 -fPIC
JULIA = julia

libkfunc.so: kfunc.o
	$(CXX) -shared -undefined dynamic_lookup $^ -o $@

kfunc.o: kfunc.cpp kfunc.h
	$(CXX) -c kfunc.cpp -o $@

clean:
	@rm -f libkfunc.so kfunc.o
