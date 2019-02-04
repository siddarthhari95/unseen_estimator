

def cartesian_product(lists):  # lists can really be any sequences
    print lists
    if any([not l for l in lists]):  # c.p. is empty if any list is empty
        return
    n = len(lists)
    indexes = [0] * n
    while True:
        print tuple(lists[i][indexes[i]] for i in xrange(n))  # currently indexed element of each list
        # update indexes
        for i in xrange(n-1, -1, -1):  # loop through indexes from back
            if indexes[i] < len(lists[i]) - 1:      # stop at first index that can be incremented ...
                indexes[i] += 1                     # ... increment it ...
                indexes[i+1:n] = [0] * (n - i - 1)  # ... reset all succeeding indexes to 0
                break
        else:  # no index could be incremented -> end
            break

import scipy.optimize as optimize

def f(params, b, x):
    # print(params)  # <-- you'll see that params is a NumPy array
    a = params # <-- for readability you may wish to assign names to the component variables
    c = a+1+b+x
    print c
    return a*2 + c

initial_guess = 0.5
b = 3
x = 100
result = optimize.minimize(f, initial_guess, args=(b, x), bounds=((17,2000),))
print result
if result.success:
    fitted_params = result.x
    print("fitted_params", fitted_params)
else:
    raise ValueError(result.message)


#  a = cartesian_product([[0,1,2,3,4],[0,1],[0,1,2,3]])
# print a
# print "cc"