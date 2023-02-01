#CC = icc
CC = gcc
CFLAGS = -O3
#CFLAGS =  -xHost -O3 
LIBS = -L./flexfloat/build -lflexfloat -lm -fopenmp

all:
	$(CC) $(CFLAGS) ./scrimp/scrimp_topK.c -o ./scrimp/scrimp_topK $(LIBS)
	$(CC) $(CFLAGS) ./scrimp/scrimp.c -o ./scrimp/scrimp $(LIBS)
	$(CC) $(CFLAGS) ./scrimp/scrimp_loc.c -o ./scrimp/scrimp_loc $(LIBS)

clean:
	rm ./scrimp/scrimp
	rm ./scrimp/scrimp_topK
	rm ./scrimp/scrimp_loc
taxi:
	./scrimp/scrimp_topK taxi.txt 100 256 0.01
taxi5:
	./scrimp/scrimp_loc taxi.txt 100 72 5 0.01
taxi25:
	./scrimp/scrimp_loc taxi.txt 100 72 25 0.01
power:
	./scrimp/scrimp power_demand.txt 200 72 1
random_similarity:
	./scrimp/scrimp random_similarity.txt 50 72 1
random_anomaly:
	./scrimp/scrimp random_anomaly.txt 50 72 1
random1:
	./scrimp/scrimp random.txt 1024 72 1 1
random:
	./scrimp/scrimp_topK random.txt 1024 256 1
random2:
	./scrimp/scrimp_loc random.txt 1024 72 2 1
ric:
	./scrimp/scrimp ric.txt 10 1 1
ric2:
	./scrimp/scrimp ric.txt 10 2 2 1
e0103:
	./scrimp/scrimp_topK e0103_n180000.txt 250 256 1
e0103_1:
	./scrimp/scrimp_topK e0103_n180000.txt 250 1 1
e0103_2:
	./scrimp/scrimp_topK e0103_n10000.txt 250 1 1
e0103_3:
	./scrimp/scrimp_topK e0103_n40000.txt 250 1 1
power-MP:
	./scrimp/scrimp_topK power-MPIII-SVF_n180000.txt 1325 256 0.1
audio-MP:
	./scrimp/scrimp_topK audio-MPIII-SVD.txt 200 256 1
