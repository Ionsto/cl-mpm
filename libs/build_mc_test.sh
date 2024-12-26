#g++ $(pkg-config --cflags eigen3) test.cpp -O3 -march=native
g++ $(pkg-config --cflags eigen3) test-mc.cpp
./a.out
