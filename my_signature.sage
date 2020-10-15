#!/usr/bin/python



# if not os.path.isfile('cable_signature.py'):
#     os.system('sage --preparse cable_signature.sage')
#     os.system('mv cable_signature.sage.py cable_signature.py')
# from cable_signature import SignatureFunction, TorusCable, SIGNATURE, SIGMA



# searching for signature > 5 + #(v_i != 0) over given knot schema
def search_for_large_signature_value(knot_formula=None, limit=None,
                                     verbose=None, print_results=None):

    if limit is None:
        limit = config.limit
    if knot_formula is None:
        knot_formula = config.knot_formula
    if verbose is None:
        verbose = config.verbose
    if print_results is None:
        print_results = config.print_results

    k_vector_size = extract_max(knot_formula) + 1
    combinations = it.combinations(range(1, limit + 1), k_vector_size)
    P = Primes()
    good_knots = []

    # iterate over q-vector
    for c in combinations:
        q = [P.unrank(i + config.start_shift) for i in c]
        if config.only_slice_candidates:
            ratio = q[3] > 4 * q[2] and q[2] > 4 * q[1] and q[1] > 4 * q[0]
            if not ratio:
                if verbose:
                    print("Ratio-condition does not hold")
                continue
        cable = TorusCable(knot_formula=knot_formula, q_vector=q)
        is_big = cable.is_signature_big_for_all_metabolizers()
        print(is_big)
        if is_big:
            good_knots.append(cable.knot_description)

    return good_knots


def extract_max(string):
    numbers = re.findall(r'\d+', string)
    numbers = map(int, numbers)
    return max(numbers)


extract_max.__doc__ = \
    """
    Return:
        maximum of absolute values of numbers from given string
    Examples:
        sage: extract_max("([1, 3], [2], [-1, -2], [-10])")
        10
        sage: extract_max("3, 55, ewewe, -42, 3300, 50")
        3300
    """

"""
This script calculates signature functions for knots (cable sums).

The script can be run as a sage script from the terminal
or used in interactive mode.

A knot (cable sum) is encoded as a list where each element (also a list)
corresponds to a cable knot, e.g. a list
[[1, 3], [2], [-1, -2], [-3]] encodes
T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7).

To calculate the number of characters for which signature function vanish use
the function eval_cable_for_null_signature as shown below.

sage: eval_cable_for_null_signature([[1, 3], [2], [-1, -2], [-3]])

T(2, 3; 2, 7) # T(2, 5) # -T(2, 3; 2, 5) # -T(2, 7)
Zero cases: 1
All cases: 1225
Zero theta combinations:
(0, 0, 0, 0)

sage:

The numbers given to the function eval_cable_for_null_signature are k-values
for each component/cable in a direct sum.

To calculate signature function for a knot and a theta value, use function
get_signature_as_function_of_theta (see help/docstring for details).

About notation:
Cables that we work with follow a schema:
    T(2, q_1; 2, q_2; 2, q_4) # -T(2, q_2; 2, q_4) #
            # T(2, q_3; 2, q_4) # -T(2, q_1; 2, q_3; 2, q_4)
In knot_formula each k[i] is related with some q_i value, where
q_i = 2*k[i] + 1.
So we can work in the following steps:
1) choose a schema/formula by changing the value of knot_formula
2) set each q_i all or choose range in which q_i should varry
3) choose vector v / theata vector.
"""
