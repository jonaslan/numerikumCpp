## **************************
## *** Numerikum, WS 2013 ***
## *** Sheet 13 Problem 3 ***
## **************************

CC		= g++
SW		= -O3 -std=c++11

# **********
main.out : main.obj odeint.obj odeint.cpp odeint.h stepperdopr5.obj stepperdopr5.cpp stepperdopr5.h vdp.obj vdp.cpp vdp.h helfer.h
	@echo "Linking 'main.out'..."
	$(CC) $(SW) odeint.obj stepperdopr5.obj vdp.obj main.obj -o main.out
	rm -f *.obj
	@echo ""

# **********
main.obj : main.cpp helfer.h
	@echo "Compiling 'main.obj'..."
	$(CC) $(SW) -c main.cpp -o main.obj
	@echo ""

# **********
odeint.obj : odeint.cpp odeint.h
	@echo "Compiling 'odeint.obj'..."
	$(CC) $(SW) -c odeint.cpp -o odeint.obj
	@echo ""

# **********
stepperdopr5.obj : stepperdopr5.cpp stepperdopr5.h
	@echo "Compiling 'stepperdopr5.obj'..."
	$(CC) $(SW) -c stepperdopr5.cpp -o stepperdopr5.obj
	@echo ""

# **********
vdp.obj : vdp.cpp vdp.h
	@echo "Compiling 'vdp.obj'..."
	$(CC) $(SW) -c vdp.cpp -o vdp.obj
	@echo ""

.PHONY: clean
clean:
	@echo "Removing compiled modules..."
	rm -f *.obj *.out
