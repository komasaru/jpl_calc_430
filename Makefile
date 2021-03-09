gcc_options = -std=c++17 -Wall -O2 --pedantic-errors

jpl_calc_430: jpl_calc_430.o jpl.o
	g++92 $(gcc_options) -o $@ $^

jpl_calc_430.o : jpl_calc_430.cpp
	g++92 $(gcc_options) -c $<

jpl.o : jpl.cpp
	g++92 $(gcc_options) -c $<

run : jpl_calc_430
	./jpl_calc_430

clean :
	rm -f ./jpl_calc_430
	rm -f ./*.o

.PHONY : run clean

