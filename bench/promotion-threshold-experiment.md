bench.pbench compare -benchmark nearestneighbors -only run --virtual_run
prun nearestneighbors.heartbeat -algorithm heartbeat -proc 40 -type array_point3d -infile _data/array_point3d_plummer_medium.bin -promotion_threshold 2,5,10,15,20,25,30,35,40,60,100,1000,10000,100000 -runs 10
pplot scatter -x promotion_threshold -y exectime --xlog --yzero


prun suffixarray.heartbeat -algorithm heartbeat -proc 40 -type string -infile _data/chr22.dna.bin -promotion_threshold 2,5,10,15,20,25,30,35,40,60,100,1000,10000,100000

prun removeduplicates.heartbeat -algorithm heartbeat -proc 40 -type array_int -infile _data/array_int_random_large.bin -promotion_threshold 2,5,10,15,20,25,30,35,40,60,100,1000,10000,100000