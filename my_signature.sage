#!/usr/bin/env python
import collections


def mod_one(n):
    """This function returns the fractional part of some number."""
    if n >= 1:
        return mod_one(n - 1)
    if n < 0:
        return mod_one(n + 1)
    return n


class av_signature_function(object):
    '''
    This simple class encodes twisted and untwisted signature functions
    of knots. Since the signature function is entirely encoded by its signature
    jump, the class stores only information about signature jumps
    in a dictionary self.data.
    The dictionary stores data of the signature jump as a key/values pair,
    where the key is the argument at which the functions jumps
    and value encodes the value of the jump. Remember that we treat
    signature functions as defined on the interval [0,1).
    '''
    def __init__(self, values=[]):
        # We will store data of signature jumps here.
        self.data = collections.defaultdict(int)
        # values contain initial data of singature jumps
        for jump_arg, jump in values:
            assert 0 <= jump_arg < 1, \
                "Signature function is defined on the interval [0, 1)."
################################### what for += ???
            self.data[jump_arg] = jump

    def value(self, arg):
        # Compute the value of the signature function at the point arg.
        # This requires summing all signature jumps that occur before arg.
        assert 0 <= arg < 1, \
            "Signature function is defined on the interval [0, 1)."
        val = 0
        for jump_arg, jump in self.data.items():
            if jump_arg < arg:
                val += 2 * jump
            elif jump_arg == arg:
                val += jump
        return val

############## what for - it is == 0
    def total_sign_jump(self):
        # Total signature jump is the sum of all jumps.
        a = sum([j[1] for j in self.to_list()])
        b = sum(self.data.values())
        # print b
        assert a == b
        return sum(self.data.values())

    def to_list(self):
        # Return signature jumps formated as a list
        return sorted(self.data.items(), key=lambda x: x[0])

    def step_function_data(self):
        # Transform the signature jump data to a format understandable
        # by the plot function.
        l = self.to_list()
        vals = ([(d[0], sum(2 * j[1] for j in l[:l.index(d)+1])) for d in l] +
                [(0, self.data[0]), (1, self.total_sign_jump())])
        return vals

    def plot(self):
        # plot the signture function
        plot_step_function(self.step_function_data())

    def tikz_plot(self, file_name):
        # Draw the graph of the signature and transform it into TiKz.
        # header of the LaTeX file
        output_file = open(file_name, "w")
        output_file.write("\\documentclass[tikz]{standalone}\n")
        output_file.write("\\usetikzlibrary{datavisualization,datavisualization.formats.functions}\n")
        output_file.write("\\begin{document}\n")
        output_file.write("\\begin{tikzpicture}\n")
        data = sorted(self.step_function_data())
        output_file.write("  \\datavisualization[scientific axes,visualize as smooth line,\n")
        output_file.write("  x axis={ticks={none,major={at={")
        output_file.write(", " + str(N(data[0][0], digits=4)) + " as \\(" + str(data[0][0]) + "\\)")
        for jump_arg, jump in data[1:]:
            output_file.write(", " + str(N(jump_arg, digits=4)) + " as \\(" + str(jump_arg) + "\\)")
            output_file.write("}}}}\n")
            output_file.write("  ]\n")
            output_file.write("data [format=function]{\n")
            output_file.write("var x : interval [0:1];\n")
            output_file.write("func y = \\value x;\n")
            output_file.write("};\n")
            # close LaTeX enviroments
            output_file.write("\\end{tikzpicture}\n")
            output_file.write("\\end{document}\n")
            output_file.close()

    def __lshift__(self, shift):
        # Shift of the signature functions correspond to the rotations.
        return self.__rshift__(-shift)

    def __rshift__(self, shift):
        new_data = []
        for jump_arg, jump in self.data.items():
            new_data.append((mod_one(jump_arg + shift), jump))
        return av_signature_function(new_data)

    def __sub__(self, other):
        # we cn perform arithmetic operations on signature functions.
        return self + other.__neg__()

    def __neg__(self):
        for jump_arg in self.data.keys():
            self.data[jump_arg] *= -1
        return self

    def __add__(self, other):
        for jump_arg, jump in other.data.items():
            self.data[jump_arg] += jump
        return self

    def __str__(self):
        return '\n'.join([str(jump_arg) + ": " + str(jump)
                          for jump_arg, jump in self.data.items()])

    def __repr__(self):
        return self.__str__()


def untw_signature(k):
    # Return the signature function of the T_{2,2k+1} torus knot.
    l = ([((2 * a + 1)/(4 * k + 2), -1) for a in range(k)] +
         [((2 * a + 1)/(4 * k + 2), 1) for a in range(k + 1, 2 * k + 1)])
    # print l
    # print type(l)
    return av_signature_function(l)
