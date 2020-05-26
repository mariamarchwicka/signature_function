This script calculates signature functions for knots (cable sums).

The script can be run as a sage script from the terminal or used in interactive mode.

A knot (cable sum) is encoded as a list where each element (also a list) corresponds to a cable knot.

To calculate the number of characters for which signature function vanish use the function eval_cable_for_thetas as shown below.

sage: eval_cable_for_thetas([[1, 3], [2], [-1, -2], [-3]])

T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)
Zero cases: 1
All cases: 1225
Zero theta combinations:
(0, 0, 0, 0)

sage:


The numbers given to the function eval_cable_for_thetas are k-values for each
component/cable in a direct sum.


To calculate signature function for a knot and a theta value, use function
get_function_of_theta_for_sum as follow:

sage: signature_function_generator = get_function_of_theta_for_sum([1, 3], [2], [-1, -2], [-3])
sage: sf = signature_function_generator(2, 1, 2, 2)
sage: print sf
0: 0
5/42: 1
1/7: 0
1/5: -1
7/30: -1
2/5: 1
3/7: 0
13/30: -1
19/42: -1
23/42: 1
            17/30: 1
4/7: 0
3/5: -1
23/30: 1
4/5: 1
6/7: 0
37/42: -1
sage:

or like below:

sage: print get_function_of_theta_for_sum([1, 3], [2], [-1, -2], [-3])(2, 1, 2, 2)
0: 0
1/7: 0
1/6: 0
1/5: -1
2/5: 1
3/7: 0
1/2: 0
4/7: 0
3/5: -1
4/5: 1
5/6: 0
6/7: 0
sage: