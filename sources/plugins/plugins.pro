TEMPLATE = subdirs

SUBDIRS  = \
				rational_fitter_cgal							\
				rational_fitter_quadprog					\
#				rational_fitter_parallel					\
				rational_fitter_eigen						\
				rational_fitter_leastsquare				\	
				rational_function_chebychev				\
				rational_fitter_matlab						\
#				rational_fitter_dca							\
				nonlinear_fitter_eigen                 \
				nonlinear_fitter_ceres						\
				nonlinear_fitter_ipopt						\
				nonlinear_fitter_nlopt						\
				nonlinear_fresnel_schlick					\
				nonlinear_function_diffuse					\
				nonlinear_function_blinn					\
				nonlinear_function_retroblinn				\
				nonlinear_function_lafortune				\
				nonlinear_function_isotropic_lafortune	\
				data_merl										\
#				data_astm

