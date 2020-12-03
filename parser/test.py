from ctypes import *
import os

so_file = "./QASM_to_Sylvan.so"
my_functions = CDLL(so_file)
os.system("./QASM_to_Sylvan ./test_qasm.txt")

my_functions.read_QASM(b"./test_qasm.txt")

