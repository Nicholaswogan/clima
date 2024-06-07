cimport OpticalProperties_pxd as op_pxd
cimport ClimaRadtranWrk_pxd as rwrk_pxd

cdef extern from *:
  struct Radtran:
    pass

cdef extern void radtran_skin_temperature_wrapper(Radtran *ptr, double *bond_albedo, double *T_skin)

cdef extern void radtran_opacities2yaml_wrapper_1(Radtran *ptr, int *out_len, void **out_cp)
cdef extern void radtran_opacities2yaml_wrapper_2(Radtran *ptr, void **out_cp, int *out_len, char* out_c)

# getters and setters
cdef extern void radtran_surface_albedo_get_size(Radtran *ptr, int *dim1)
cdef extern void radtran_surface_albedo_get(Radtran *ptr, int *dim1, double *arr)
cdef extern void radtran_surface_albedo_set(Radtran *ptr, int *dim1, double *arr)

cdef extern void radtran_surface_emissivity_get_size(Radtran *ptr, int *dim1)
cdef extern void radtran_surface_emissivity_get(Radtran *ptr, int *dim1, double *arr)
cdef extern void radtran_surface_emissivity_set(Radtran *ptr, int *dim1, double *arr)

cdef extern void radtran_ir_get(Radtran *ptr, op_pxd.OpticalProperties **ptr1);

cdef extern void radtran_sol_get(Radtran *ptr, op_pxd.OpticalProperties **ptr1);

cdef extern void radtran_wrk_ir_get(Radtran *ptr, rwrk_pxd.ClimaRadtranWrk **ptr1);

cdef extern void radtran_wrk_sol_get(Radtran *ptr, rwrk_pxd.ClimaRadtranWrk **ptr1);