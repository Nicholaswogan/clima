
cdef extern void radtran_skin_temperature_wrapper(void *ptr, double *bond_albedo, double *T_skin)

cdef extern void radtran_opacities2yaml_wrapper_1(void *ptr, int *out_len, void *out_cp)
cdef extern void radtran_opacities2yaml_wrapper_2(void *ptr, void *out_cp, int *out_len, char* out_c)

# getters and setters
cdef extern void radtran_surface_albedo_get(void *ptr, double *val)
cdef extern void radtran_surface_albedo_set(void *ptr, double *val)

cdef extern void radtran_ir_get(void *ptr, void *ptr1);

cdef extern void radtran_sol_get(void *ptr, void *ptr1);

cdef extern void radtran_wrk_ir_get(void *ptr, void *ptr1);

cdef extern void radtran_wrk_sol_get(void *ptr, void *ptr1);