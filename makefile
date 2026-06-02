CFLAGS      = -Wall -O3 $(EXTRA_CFLAGS)
LFLAGS      = -Wall -O3
CC          = gcc

OBJ         = obj/heap.o obj/boastar.o obj/graph.o obj/main_boa.o obj/pool.o
OBJ_BOD     = obj/heap.o obj/bod.o obj/graph.o obj/main_bod.o
OBJ_NAMOADR = obj/heap.o obj/namoadr.o obj/graph.o obj/main_namoadr.o
OBJ_BCP_BOA = obj/heap.o obj/boastar.o obj/bcp_boastar.o obj/graph.o obj/main_bc_boastar.o obj/pool.o
OBJ_BENCH   = obj/heap.o obj/boastar.o obj/bcp_boastar.o obj/graph.o obj/benchmark.o obj/pool_bench.o

all: boa bod namoadr bcp_boastar benchmark

boa: $(OBJ)
	$(CC) $(LFLAGS) -o boa $(OBJ) -lm

bod: $(OBJ_BOD)
	$(CC) $(LFLAGS) -o bod $(OBJ_BOD) -lm

namoadr: $(OBJ_NAMOADR)
	$(CC) $(LFLAGS) -o namoadr $(OBJ_NAMOADR) -lm

bcp_boastar: $(OBJ_BCP_BOA)
	$(CC) $(LFLAGS) -o bcp_boastar $(OBJ_BCP_BOA) -lm

benchmark: $(OBJ_BENCH)
	$(CC) $(LFLAGS) -o benchmark $(OBJ_BENCH) -lm

obj/%.o: src/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f obj/*.o boa bod namoadr bcp_boastar benchmark

rebuild: clean all

test-BAY:
	./benchmark Maps/BAY-road-d.txt 10 30