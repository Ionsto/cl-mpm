#g++ $(pkg-config --cflags eigen3) test.cpp -O3 -march=native
#$CC -fPIC -shared $(pkg-config --cflags eigen3) test.cpp -O3 -march=native -o kirchoff.so
#$CXX -fPIC -I ~/eigen/ test.cpp -O3 -march=native -o kirchoff.so

#$CXX -fPIC -I ~/eigen/ test.cpp -O3 -march=native -shared -o kirchoff.so
export CXX=g++
$CXX -fPIC $(pkg-config --cflags eigen3) test.cpp -O3 -march=native -shared -o kirchoff.so
#$CXX -fPIC $(pkg-config --cflags eigen3) test.cpp -O3 -march=native -shared -o kirchoff.so
