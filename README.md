## MultiParty Communications on Ring

### ***environment***
g++:  sudo apt-get install build-essential  
m4:  sudo apt-get install m4  
gmp:  cd gmp's catagory; ./configure; make; make install; sudo make install(ls /usr/local/include/; ls /usr/local/lib/)  
ntl:  cd ntl's catagory; sudo ldconfig; ./configure; make; make install; sudo make install(ls /usr/local/include/; ls /usr/local/lib/)  

  
### ***run command***
g++ -g -O2 -std=c++11 -pthread -march=native correspond.cpp -o out -lntl -lgmp -lm
