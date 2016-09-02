ifdef opt
	copt = -O$(opt)
else
	copt =
endif
#clean for testing makefile
main : clean main.obj xorshift.obj moments.obj
	g++ -Wall $(copt) main.obj xorshift.obj moments.obj -o main

main.obj : main.cpp
	g++ -c -Wall $(copt) main.cpp -o main.obj

xorshift.obj : xorshift.h xorshift.cpp
	g++ -c -Wall $(copt) xorshift.cpp -o xorshift.obj

moments.obj : moments.h moments.cpp
	g++ -c -Wall $(copt) moments.cpp -o moments.obj

.PHONY : clean
clean :
	rm -f *.obj
