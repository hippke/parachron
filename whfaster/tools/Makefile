REB=/home/claude/rebound/src
CFLAGS=-O3 -march=native -I$(REB) -Wall
LDFLAGS=-L$(REB) -lrebound -Wl,-rpath,$(REB) -lm
%: %.c ss_ics.h
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
