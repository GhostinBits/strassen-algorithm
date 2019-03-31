# strassen-algorithm
## with a strange addressing way
Implementation of the algorithm as described in Introduction to Algorithm with C++ (required by my instructor, though it's not my favourite, reasons given below)

2D pointer as parametre is not supported by C++, where I struggle to find the reasons since one could use 1D pointer as a parametre without providing the array's length and the same technique couldn't be used with 2D pointer. Providing the array's length is no where near elegence. Why should people and how could users compile the code each time they want to compute some matrices? And C++ doesn't support VLA as well.

That's why the implemention doesn't utilise 2D arrays to reprensent matrices. Instead, the implemention avoids providing the length by using 1D pointer and simply uses algebra to address.
