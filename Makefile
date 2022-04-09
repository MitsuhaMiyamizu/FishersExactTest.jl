CXX = clang++
CXXFLAGS = -std=c++17 -fPIC
JULIA = julia

test: libkfunc.so test.jl
	$(JULIA) test.jl

libkfunc.so: kfunc.o
	$(CXX) -shared -undefined dynamic_lookup $^ -o $@

kfunc.o: kfunc.cpp kfunc2.h
	$(CXX) -c kfunc.cpp -o $@

clean:
	@rm -f libkfunc.so kfunc.o
