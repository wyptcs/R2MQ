This code creates an encoding data structure to compute R2MQ(i,j) on an array of unsigned integers.
Contact : wypark2510@gmail.com

**data structures for answering RT2Q**
<br/>
<br/>
Description:<br/>
This is an RT2Q compressed data structure on 'RT2Qquery' folder. This experiment calculates query time and space usage. We used Range Min-Max Tree of Sadakane and Navarro[[1] on SDSL library[2]. In order to reduce the size, our Range Min-Max Tree uses minimum field only. Also, this includes three RMQ structures for comparison, based on the method of Ferrada and Navarro[3], Fischer and Heun[4], BGHL[5]
We included our test example which creates experimental results of our paper, and toy example based on Figure 1 of our paper.

Make:<br/>

For space usage and query time test<br/>

To make the library give the command make and this will create the lib: 'DFUDSR2MQ.a'. To test our toy example, give the command 'make tutorial', and to test our experiment, give the command 'make example'(it is for ours and our implementation based on Fischer and Heun[4].<br/>
To test other comparison test cases, give the command 'make dccdfuds' 'make dccbp' 'make sea'<br/>
summary : make && make example (ours)<br/>
make && make sea [2]<br/>
make && make dccdfuds [3]<br/>
make && make dccbp [3]<br/>

Compile:<br/>
To use the library you must complie your program linking 'DFUDSR2MQ.a' and include the file for what you are using.<br/>
For example to compile the file example.cpp (included here) we will run:<br/>
g++ -std=c++11 -O3 -DNDEBUG example.cpp -o example DFUDSR2MQ.a -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64<br/>
This binary have to receive 4 parameter.<br/>
1.- mode : execute which Top-2 structure, 0 : our encoding, 3 : Fischer, Heun Indexing<br/>
2.- path : integer array file path, integer should be splited using '\n'<br/>
3.- query size : gives query size j-i+1 on query range [i,j]<br/>
4.- using space(our encoding only) : it gives depth/ldepth structure size. For example, parameter 0.25 means we use 0.25n bits of depth/ldepth structure.<br/>
<br/>
<br/>
For example, when we use 0.25n bits for depth/ldepth structure size for test query size 1000 and use our encoding for "RANDOM5.txt" file:
./example 0 RANDOM5.txt 1000 0.25<br/>
To test [4] at the same file and query size, use ./example 3 RANDOM5.txt 1000<br/>
Also, to test [2] at the same file and query size, use  ./SEATEST RANDOM5.txt 1000 command.<br/>
to test [3], use ./dcctest RANDOM5.txt 1000 on dfuds mode or ./dccbptest RANDOM5.txt 1000 on bp mode<br/>
<br/>
<br/>
<br/>

**space-efficient encodings for RT2Q**
<br/>
<br/>
Description:<br/>
This is an RT2Q compressed data structure on 'spaceefficient' folder. This experiment calculates encoding time and space usage. We used our data structure with DAG structure[6] and Optimal structure[7].

Make:<br/>
For encoding time and encoding space test<br/>
Just give the command make and it will make binary 'spacetest'<br/>

Compile:<br/>
To conduct experiment, 'spacetest' takes three parameters.<br/>
1.- path : integer array file path, integer should be splited using '\n'<br/>
2.- optmode : when constructing DAG, to use optimization scheme give 1 else give 0.<br/>
3.- revmode : when constructing our data structure, value 1 means reverse the original array else value 0.<br/>
For example, to test using optimizing and not using revmode with RANDOM5.txt, use ./spacetest RANDOM5.txt 1 0.<br/>
<br/>
<br/>

**References**<br/>
[1]. K. Sadakane and G. Navarro. Fully-Functional Static and Dynamic Succinct Trees. ACM Transactions on Algorithms 10(3):article 16, 2014.<br/>
[2]. Simon Gog, Timo Beller, Alistair Moffat, and Matthias Petri, “From theory to practice:Plug  and  play  with  succinct  data  structures,”   in13th International Symposium onExperimental Algorithms, (SEA 2014), 2014, pp. 326–337.<br/>
[3]. Hector Ferrada and Gonzalo Navarro, “Improved range minimum queries,”J. DiscreteAlgorithms, vol. 43, pp. 72–80, 2017.<br/>
[4]. J. Fischer and V. Heun. Space-efficient preprocessing schemes for range minimum queries on static arrays. SIAM Journal on Computing, 40(2):465–492, 2011.<br/>
[5]. Baumstark N, Gog S, Heuer T, Labeit J. Practical range minimum queries revisited. In: Proceedings of the 16th International Symposium on Experimental Algorithms (SEA 2017); 2017; London. UK.<br/>
[6]. S. Jo, R. Lingala, S. R. Satti, Encoding two-dimensional range top-k queries, in: CPM 2016, 2016, pp. 3:1–3:11.
[7]. P. Gawrychowski, P. K. Nicholson, Optimal encodings for range top-k, selection, and min-max, in: ICALP 2015,
Proceedings, Part I, 2015, pp. 593–604.
