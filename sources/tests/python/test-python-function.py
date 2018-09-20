import alta
import sys

# Test a nonlinear_function
# Open a function plugin, then save the function and load it again
print('Loading a nonlinear_function plugin and testing saving/loading')
f1 = alta.get_function('nonlinear_function_diffuse',
                       alta.parameters(3, 3,
                                       alta.input_parametrization.RUSIN_TH_PH_TD,
                                       alta.output_parametrization.RGB_COLOR))
print("f1 = ", f1)
f1.save('test-diffuse.func')

# Test a rational_function 
# Open a function plugin, then save the function and load it again
print('Loading a rational_function plugin and testing saving/loading')
f2 = alta.get_function('rational_function_legendre',
                       alta.parameters(3, 3,
                                       alta.input_parametrization.RUSIN_TH_PH_TD,
                                       alta.output_parametrization.RGB_COLOR))
print("f2 = ", f2)
f2.save('test-rat.func')

sys.exit(0)
