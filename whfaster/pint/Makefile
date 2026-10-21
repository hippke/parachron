REB=/home/claude/rebound/src
CFLAGS=-O3 -march=native -I$(REB) -Wall
LDFLAGS=-L$(REB) -lrebound -Wl,-rpath,$(REB) -lm -lpthread
pint_test: pint_test.c pint_coarse.c pint_coarse.h
	$(CC) $(CFLAGS) pint_test.c pint_coarse.c -o $@ $(LDFLAGS)
pint_run: pint_run.c pint_coarse.c pint_coarse.h
	$(CC) $(CFLAGS) pint_run.c pint_coarse.c -o $@ $(LDFLAGS)
libpint.so: pint_coarse.c pint_coarse.h
	$(CC) $(CFLAGS) -shared -fPIC pint_coarse.c -o $@ $(LDFLAGS)
pint_long: pint_long.c pint_coarse.c pint_coarse.h
	$(CC) $(CFLAGS) pint_long.c pint_coarse.c -o $@ $(LDFLAGS)
pint_pipe: pint_pipe.c pint_coarse.c pint_coarse.h
	$(CC) $(CFLAGS) pint_pipe.c pint_coarse.c -o $@ $(LDFLAGS)
ser_pert: ser_pert.c pint_coarse.c pint_coarse.h
	$(CC) $(CFLAGS) ser_pert.c pint_coarse.c -o $@ $(LDFLAGS)
