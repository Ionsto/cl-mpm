#g++ $(pkg-config --cflags eigen3) test.cpp -O3 -march=native
g++ -fPIC -shared $(pkg-config --cflags eigen3) test.cpp -O3 -march=native -o kirchoff.so
#./a.out
