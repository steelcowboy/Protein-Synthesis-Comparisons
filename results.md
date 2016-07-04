I ran each program 5 times and took the average using:

time <program> -ni > /dev/null

Note that C++_m means multithreaded, using 4 threads.
Additionally, all cout statements were commented out
in this version to avoid slowdown. For comparison purposes,
the standard version also was tested with cout statements
commented out. The asterisk means it prints to cout.

Results
--------------
Python: 1.8582
C++*:   0.6516
C++:    0.5088
C++_m:  0.2554



On a Raspberry Pi 2 the results are more distinctive:

Results: 
-------------
Python: 23.3114
C++*:   5.0008
C++:    3.5061
C++_m:  1.0146


As of Version 0.5 the Python version functions exactly as
the C++ version, using a library with all of the same function
calls. I also ran the Cython compiler on just the module with
pure Python to compare the results.

Results will be added in V0.6 with synthesis.py cythonized and
testing the C++ interface with the object built by Cython. If
the impact on speed is negligible, I may consider using Cython
as the primary language for the library


Results: 
-------------
Python with normal module: 1.3534
Python with cython module: 0.8616
~~Python with cythonized module:~~
~~C++ compiled with cythonized module:~~


