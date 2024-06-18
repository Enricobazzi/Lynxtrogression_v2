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

def p0_6_2():
    p0 = [0.18370924973853417, 0.003457059809692882, 0.29853228460362985, 0.1628426520848679, 0.030399222195332302, 0.011001828266108864, 1.071051337716487, 2.867613590978363, 0.0034764198065980905, 0.01855232203947269, 0.017514079791411515, 0.0012357341515685301, 0.0007780542768669069, 0.0016140472420234388, 0.01078524915659998, 0.005544824954306601, 0.983202690396147, 0.9019266299149103]
    return p0

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

def p0_12_9():
    p0 = [0.0034143275463893994, 0.017342790687784326, 0.1675449264510798, 0.35093849895789375, 0.18687814605555078, 0.12387824327714872, 0.0, 0.10331917142947933, 0.0015918914318867083, 0.08823403167785245, 0.012110421222315683, 1.1635631820989543, 0.13059537751982772, 0.003571521635780008, 0.09042110097164137, 0.0431890898431185, 0.004258757033436549, 0.133863703771947]
    return p0

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

def p0_20_7():
    p0 = [0.05776609464887263, 0.001, 0.01714444072423897, 0.1316895716275381, 0.012360070259564848, 0.006737484793477517, 1.8431197984611976, 3.0596060259498925, 0.12931945297267977, 0.0796289079962728, 0.026070392553972605, 0.7825871436408913, 1.4208734865261692, 0.00655424769194445, 0.020062330959718308, 0.006699737738753766, 0.5034922972251962, 0.5454925500503207]
    return p0

def pnames_20_7():
	pnames = ['t1', 'nu11', 'nu11_1', 't2', 'nu21', 'nu22', 'm2_12', 'm2_21', 't3', 'nu31', 'nu32', 'm3_12', 'm3_21', 't4', 'nu41', 'nu42', 'm4_12', 'm4_21']
	return pnames
