#ifdef __cplusplus
 extern "C"
{
#endif

typedef struct gate_info_s {
    char* gate;
    double angle;
    uint32_t target;
    uint32_t controls[3];
} gate_info_t;

gate_info_t* parse_qasm_file(char *filepath);

#ifdef __cplusplus
}
#endif