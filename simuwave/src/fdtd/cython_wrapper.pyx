cdef extern from "fdtd3d.h":
	object run_fdtd(object ipt)


def run_fdtd_python(ipt):
	"""Executes the C code containing the actual FDTD algorithm. """
	return run_fdtd(ipt)
