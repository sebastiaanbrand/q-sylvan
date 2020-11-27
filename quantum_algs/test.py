from ctypes import *
import os

so_file = "./QASM_to_Sylvan.so"
my_functions = CDLL(so_file)

my_functions.read_QASM(b"./test_qasm.txt")

