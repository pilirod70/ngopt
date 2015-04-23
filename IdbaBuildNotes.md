On an ubuntu 10.10 system with gcc 4.4.5, the following will build idba:
```
wget http://hku-idba.googlecode.com/files/idba-0.18.tar.gz
tar xvzf idba-0.18.tar.gz
cd idba-0.18
./configure --prefix=$HOME
```

Now one must edit the file src/tools/spliceGraph.cpp
Line 287 contains a typo, where a variable called "index" is used but it should instead be "i".  So make that change then:

```
make install
```