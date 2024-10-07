import dadi

def model_func_6_2(params, ns, pts):
	t1, nu11, nu11_2, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu11)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu2_func = lambda t: nu11_2 * (nu22 / nu11_2) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu21, nu2=nu2_func, m12=m2_12, m21=m2_21)
	nu1_func = lambda t: nu21 * (nu31 / nu21) ** (t / t3)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu1_func, nu2=nu32, m12=m3_12, m21=m3_21)
	nu1_func = lambda t: nu31 * (nu41 / nu31) ** (t / t4)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu1_func, nu2=nu42, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_6_2():
	pnames = ['t1', 'nu11', 'nu11_2', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames

def model_func_12_9(params, ns, pts):
	t1, nu11, nu11_1, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu11)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu11_1 * (nu21 / nu11_1) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu1_func, nu2=nu22, m12=m2_12, m21=m2_21)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu31, nu2=nu32, m12=m3_12, m21=m3_21)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu41, nu2=nu42, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_12_9():
	pnames = ['t1', 'nu11', 'nu11_1', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames

def model_func_20_7(params, ns, pts):
	t1, nu11, nu11_1, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	nu1_func = lambda t: _Nanc_size * (nu11 / _Nanc_size) ** (t / t1)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu1_func)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu11_1 * (nu21 / nu11_1) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu1_func, nu2=nu22, m12=m2_12, m21=m2_21)
	nu1_func = lambda t: nu21 * (nu31 / nu21) ** (t / t3)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu1_func, nu2=nu32, m12=m3_12, m21=m3_21)
	nu2_func = lambda t: nu32 * (nu42 / nu32) ** (t / t4)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu41, nu2=nu2_func, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_20_7():
	pnames = ['t1', 'nu11', 'nu11_1', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames

def model_func_34_7(params, ns, pts):
	t1, nu11, nu11_1, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu11)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu11_1 * (nu21 / nu11_1) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu1_func, nu2=nu22, m12=m2_12, m21=m2_21)
	nu2_func = lambda t: nu22 * (nu32 / nu22) ** (t / t3)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu31, nu2=nu2_func, m12=m3_12, m21=m3_21)
	nu2_func = lambda t: nu32 * (nu42 / nu32) ** (t / t4)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu41, nu2=nu2_func, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_34_7():
	pnames = ['t1', 'nu11', 'nu11_1', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames

def model_func_38_4(params, ns, pts):
	t1, nu11, nu11_1, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu11)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu11_1 * (nu21 / nu11_1) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu1_func, nu2=nu22, m12=m2_12, m21=m2_21)
	nu1_func = lambda t: nu21 * (nu31 / nu21) ** (t / t3)
	nu2_func = lambda t: nu22 * (nu32 / nu22) ** (t / t3)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu1_func, nu2=nu2_func, m12=m3_12, m21=m3_21)
	nu1_func = lambda t: nu31 * (nu41 / nu31) ** (t / t4)
	nu2_func = lambda t: nu32 * (nu42 / nu32) ** (t / t4)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu1_func, nu2=nu2_func, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_38_4():
	pnames = ['t1', 'nu11', 'nu11_1', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames

def model_func_30_1(params, ns, pts):
	t1, nu11, nu11_2, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu11)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu2_func = lambda t: nu11_2 * (nu22 / nu11_2) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu21, nu2=nu2_func, m12=m2_12, m21=m2_21)
	nu1_func = lambda t: nu21 * (nu31 / nu21) ** (t / t3)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu1_func, nu2=nu32, m12=m3_12, m21=m3_21)
	nu2_func = lambda t: nu32 * (nu42 / nu32) ** (t / t4)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu41, nu2=nu2_func, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_30_1():
	pnames = ['t1', 'nu11', 'nu11_2', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames

def model_func_18_7(params, ns, pts):
	t1, nu11, nu11_1, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu11)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu11_1 * (nu21 / nu11_1) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu1_func, nu2=nu22, m12=m2_12, m21=m2_21)
	nu2_func = lambda t: nu22 * (nu32 / nu22) ** (t / t3)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu31, nu2=nu2_func, m12=m3_12, m21=m3_21)
	nu1_func = lambda t: nu31 * (nu41 / nu31) ** (t / t4)
	nu2_func = lambda t: nu32 * (nu42 / nu32) ** (t / t4)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu1_func, nu2=nu2_func, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_18_7():
	pnames = ['t1', 'nu11', 'nu11_1', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames

def model_func_18_10(params, ns, pts):
	t1, nu11, nu11_2, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu11)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu2_func = lambda t: nu11_2 * (nu22 / nu11_2) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu21, nu2=nu2_func, m12=m2_12, m21=m2_21)
	nu2_func = lambda t: nu22 * (nu32 / nu22) ** (t / t3)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu31, nu2=nu2_func, m12=m3_12, m21=m3_21)
	nu1_func = lambda t: nu31 * (nu41 / nu31) ** (t / t4)
	nu2_func = lambda t: nu32 * (nu42 / nu32) ** (t / t4)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu1_func, nu2=nu2_func, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_18_10():
	pnames = ['t1', 'nu11', 'nu11_2', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames

def model_func_12_6(params, ns, pts):
	t1, nu11, nu11_1, nu11_2, t2, nu21, nu22, m2_12, m2_21, t3, nu31, nu32, m3_12, m3_21, t4, nu41, nu42, m4_12, m4_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, T=t1, nu=nu11)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu11_1 * (nu21 / nu11_1) ** (t / t2)
	nu2_func = lambda t: nu11_2 * (nu22 / nu11_2) ** (t / t2)
	phi = dadi.Integration.two_pops(phi, xx, T=t2, nu1=nu1_func, nu2=nu2_func, m12=m2_12, m21=m2_21)
	nu1_func = lambda t: nu21 * (nu31 / nu21) ** (t / t3)
	phi = dadi.Integration.two_pops(phi, xx, T=t3, nu1=nu1_func, nu2=nu32, m12=m3_12, m21=m3_21)
	nu2_func = lambda t: nu32 * (nu42 / nu32) ** (t / t4)
	phi = dadi.Integration.two_pops(phi, xx, T=t4, nu1=nu41, nu2=nu2_func, m12=m4_12, m21=m4_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

def pnames_12_6():
	pnames = ['t1', 'nu11', 'nu11_1', 'nu11_2', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames
