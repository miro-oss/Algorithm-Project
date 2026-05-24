CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall

all: generator resequencer

generator: generator.cpp
	$(CXX) $(CXXFLAGS) -o generator generator.cpp

resequencer: resequencer.cpp
	$(CXX) $(CXXFLAGS) -o resequencer resequencer.cpp

run: all
	mkdir -p data
	./generator
	./resequencer

clean:
	rm -f generator resequencer
	rm -rf data/
