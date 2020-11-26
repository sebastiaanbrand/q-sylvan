        	# hello
from ctypes import *

so_file = "./QASM_to_Sylvan.so"    # what to do   
my_functions = CDLL(so_file)

my_functions.read_QASM(b"./test_qasm.txt")

