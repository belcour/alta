import alta

fail = False

args = alta.arguments({'hello' : 'world'})
if args['hello'] != 'world':
    fail = True

args = alta.arguments({'hello' : 'world',
                       'foo'   : 'bar',
                       'param' : 'STARK_2D'})
if args['param'] != 'STARK_2D':
    fail = True

if fail:
    print("Testing python 'arguments' interface failed!")
    exit(1)
else:
    print("Testing python 'arguments' interface passed!")
    exit(0)
