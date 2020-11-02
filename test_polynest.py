# The naming conventions are important so we can EASILY pass things back and forth.
# You name your stuff and put the inputs in order so they match these test files.
# For each week (chapter) you will create one file with all implementations.
# You will name that file numerics#_lastname.py  where
#    you replace # with the chapter number and
#    you replace lastname with your lastname.'
# For this first day, we have our own file for chapter 0 with only one function in it.
# That function should be called like this:
#        output = polynest(x, coeffvec, basepointvec)
#
#
# Change the next line to import your own numerics0 file:
import numerics0_he as num
import numpy as np

coeff = [1, 2, 1, 3]
bspts = [2, 2, 2]

out = num.polynest(1,coeff)

print(out)
print('The value should be 7')

outbs = num.polynest(1,coeff,bspts)
print(outbs)
print('The value should be -3')

# You don't need to do these things if you don't want to do them.
badinput = np.array([coeff,coeff])
print(badinput.ndim)
#
outbad = num.polynest(1,badinput)

# I encourage to make up some polynomials and evaluate them.





























