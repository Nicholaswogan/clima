
cdef extern void radtran_skin_temperature_wrapper(void *ptr, double *bond_albedo, double *T_skin)

# getters and setters
cdef extern void radtran_surface_albedo_get(void *ptr, double *val)
cdef extern void radtran_surface_albedo_set(void *ptr, double *val)

cdef extern void radtran_ir_get(void *ptr, void *ptr1);

cdef extern void radtran_sol_get(void *ptr, void *ptr1);

cdef extern void radtran_wrk_ir_get(void *ptr, void *ptr1);

cdef extern void radtran_wrk_sol_get(void *ptr, void *ptr1);