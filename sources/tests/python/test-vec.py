import alta
import numpy

fail = 0

# Init the vector using an input array
# Also test if the _getitem_ function works
x = alta.vec([1.0, 1.0, 1.0])
if x[0] != 1.0:
    fail += 1

# Check if _setitem_ function works
x[0] = 0.0
if x[0] != 0.0:
    fail += 1

# Vector arithmetic
y = x + x
z = x - x
if y[0] != 2.0 and z[0] != 0.0:
    fail += 1

# Testing str(vec)
print "Testing 'str' method on variable 'y': " + str(y)

# Testing the length of the vector
if len(y) != 3:
    fail += 1

# Conversion to a list
l = list(y)
if len(l) != len(y):
    fail += 1
for i in range(0,len(y)):
    if l[i] != y[i]:
        fail += 1

# Conversion to a numpy array
na = numpy.array(y)
for i in range(0,len(y)):
    if na[i] != y[i]:
        fail += 1


# Print the test result
if fail > 0:
    print "Testing python 'vec' interface failed!"
    exit(1)
else:
    print "Testing python 'vec' interface passed!"
    exit(0)
