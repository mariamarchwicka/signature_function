#!/usr/bin/env python

import collections
import sys


def mod_one(n):
    """This function returns the fractional part of some number."""
    n -= int(n)
    if n < 0:
        n += 1
    return n


class av_signature_function(object):
    """
    This simple class encodes twisted and untwisted signature functions
    of knots. Since the signature function is entirely encoded by its signature
    jump, the class stores only information about signature jumps
    in a dictionary self.data.
    The dictionary stores data of the signature jump as a key/values pair,
    where the key is the argument at which the functions jumps
    and value encodes the value of the jump. Remember that we treat
    signature functions as defined on the interval [0,1).
    """
    def __init__(self, values=[]):
        # We will store data of signature jumps here.
        self.data = collections.defaultdict(int)
        # values contain initial data of singature jumps
        for jump_arg, jump in values:
            assert 0 <= jump_arg < 1, \
                "Signature function is defined on the interval [0, 1)."
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

    def sum_of_values(self):
        # Total signature jump is the sum of all jumps.
        a = sum([j[1] for j in self.to_list()])
        b = sum(self.data.values())
        # print b
        assert a == b
        assert a == 0
        return sum(self.data.values())

    def sum_of_absolute_values(self):
        return sum([abs(i) for i in self.data.values()])

    def double_cover(self):
        new_data = []
        for jump_arg, jump in self.data.items():
            new_data.append((mod_one(jump_arg/2), jump))
            new_data.append((mod_one(1/2 + jump_arg/2), jump))
        return av_signature_function(new_data)


    def to_list(self):
        # Return signature jumps formated as a list
        return sorted(self.data.items(), key=lambda x: x[0])

    def step_function_data(self):
        # Transform the signature jump data to a format understandable
        # by the plot function.
        l = self.to_list()
        vals = ([(d[0], sum(2 * j[1] for j in l[:l.index(d)+1])) for d in l] +
                [(0, self.data[0]), (1, self.sum_of_values())])
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
        new_one = av_signature_function()
        new_data = collections.defaultdict(int)
        for jump_arg, jump in other.data.items():
            new_data[jump_arg] = jump + self.data.get(jump_arg, 0)
            try:
                int(jump_arg)
            except:
                print jump_arg
        for jump_arg, jump in self.data.items():
            if jump_arg not in new_data.keys():
                new_data[jump_arg] = self.data[jump_arg]

        new_one.data = new_data
        return new_one

    def __str__(self):
        return '\n'.join([str(jump_arg) + ": " + str(jump)
                          for jump_arg, jump in sorted(self.data.items())])

    def __repr__(self):
        return self.__str__()

# 9.8
# ksi = exp( (2 PI * i) / (2k + 1))
# blanchfield = lambda_even + lambda_odd

def get_twisted_signature_function(k_n, theta):
    results = []
    k = abs(k_n)

    ksi = 1/(2 * k + 1)
    # lambda_odd (theta + e) % 2 == 0:
    for e in range(1, k + 1):
        if (theta + e) % 2 != 0:
            results.append((e * ksi, 1 * sgn(k_n)))
            results.append((1 - e * ksi, -1 * sgn(k_n)))
    # lambda_even
    # print "normal"
    for e in range(1, theta):
        if (theta + e) % 2 == 0:
            results.append((e * ksi, 1 * sgn(k_n)))
            results.append((1 - e * ksi, -1 * sgn(k_n)))
    # print "reversed"
    for e in range(theta + 1, k + 1):
        if (theta + e) % 2 != 0:
            continue
        results.append((e * ksi, -1 * sgn(k_n)))
        results.append((1 - e * ksi, 1 * sgn(k_n)))
    return av_signature_function(results)

def get_blanchfield(t, k):
    p = 2
    q = 2 * k + 1
    sigma_set = get_sigma_set(p, q)
    sigma = len(sigma_set) - 2 * len([z for z in sigma_set if t < z < 1 + t])
    return sigma

def get_sigma_set(p, q):
    sigma_set = set()
    for i in range(1, p):
        for j in range(1, q):
            sigma_set.add(j/q + i/p)
    return sigma_set

# Bl_theta(K'_(2, d) = Bl_theta(T_2, d) + Bl(K')(ksi_l^(-theta) * t) + Bl(K')(ksi_l^theta * t)

def get_cable_signature_as_theta_function(*arg):
    def signture_function(theta):
        if theta > abs(arg[-1]):
            print "k for pattern is " + str(arg[-1])
            print "theta shouldn't be larger than this"
            return None
        if theta == 0:
            cable_signature = get_untwisted_signutere_function(arg[-1])
        else:
            cable_signature = get_twisted_signature_function(arg[-1], theta)

        for i, k_i in enumerate(arg[:-1][::-1]):
            k = abs(k_i)
            ksi = 1/(2 * k + 1)
            power = 2^i
            a = get_untwisted_signutere_function(k_i)
            shift = theta * ksi * power
            b = a >> shift
            c = a << shift
            for _ in range(i):
                b = b.double_cover()
                c = c.double_cover()
            b += c
            cable_signature += b
        return cable_signature
    return signture_function

def get_untwisted_signutere_function(*arg):
    signture_function = av_signature_function([(0, 0)])
    for k_i in arg:
        k = abs(k_i)
        # Return the signature function of the T_{2,2k+1} torus knot.
        l = ([((2 * a + 1)/(4 * k + 2), -1 * sgn(k_i)) for a in range(k)] +
             [((2 * a + 1)/(4 * k + 2), 1 * sgn(k_i)) for a in range(k + 1, 2 * k + 1)])
        signture_function += av_signature_function(l)
    return signture_function

def get_function_of_theta_for_sum(*arg):
    def signture_function_for_sum(*thetas):
        if len(thetas) != len(arg) - 1:
            print "For each cable one theta value should be given"
            return None
        signature_function = get_untwisted_signutere_function(*arg[0])
        for i, knot in enumerate(arg[1:]):
            signature_function += (get_cable_signature_as_theta_function(*knot))(thetas[i])
        return signature_function
    return signture_function_for_sum

def first_sum(k_0, k_1, k_2, k_3):
    F = get_function_of_theta_for_sum([k_3, -k_2], [-k_0, -k_1, -k_3], [k_0, k_1, k_2])
    for theta_0 in range(k_3 + 1):
        for theta_1 in range(k_2 + 1):
            f = F(theta_0, theta_1)
            f.sum_of_values()
            if f.sum_of_absolute_values() != 0 and theta_1 + theta_0 == 0:
                    print 4 * "\n"
                    print "OJOJOJOJJOOJJOJJ!!!!!!!!!!"
                    print k_0, k_1, k_2, k_3
                    print theta_0, theta_1

            if f.sum_of_absolute_values() == 0 and theta_1 + theta_0 != 0:
                # print "HURA"
                # print k_0, k_1, k_2, k_3
                # print theta_0, theta_1
                if k_2 != k_3 or theta_0 != theta_1:
                    print 4 * "\n"
                    print " SUPER!!!!!!!!!!"
                    print k_0, k_1, k_2, k_3
                    print theta_0, theta_1

def second_sum(k_0, k_1, k_2, k_3, k_4):
    F = get_function_of_theta_for_sum([], [k_0, k_1, k_2], [k_3, k_4], [-k_0, -k_3, -k_4], [-k_1, -k_2])
    for theta_0 in range(k_2 + 1):
        for theta_1 in range(k_4 + 1):
            for theta_2 in range(k_4 + 1):
                for theta_3 in range(k_2 + 1):
                    f = F(theta_0, theta_1, theta_2, theta_3)
                    if f.sum_of_absolute_values() != 0 and theta_1 + theta_0 + theta_3 + theta_2 == 0:
                            print 4 * "\n"
                            print "2 OJOJOJOJJOOJJOJJ!!!!!!!!!!"
                            print k_0, k_1, k_2, k_3, k_4
                            print theta_0, theta_1, theta_2, theta_3

                    if f.sum_of_absolute_values() == 0 and theta_1 + theta_0 + theta_3 + theta_2 != 0:
                        # print "HURA"
                        # print k_0, k_1, k_2, k_3
                        # print theta_0, theta_1
                        if k_2 != k_3 or theta_0 != theta_1:
                            print 4 * "\n"
                            print "2 SUPER!!!!!!!!!!"
                            print k_0, k_1, k_2, k_3, k_4
                            print theta_0, theta_1, theta_2, theta_3

def third_sum(k_0, k_1, k_2, k_3, k_4, k_5, k_6, k_7, k_8):
    F = get_function_of_theta_for_sum([], [k_0, k_1, k_2], [k_3, k_4], [-k_5, -k_6, -k_7], [-k_8, -k_8])
    for theta_0 in range(k_2 + 1):
        for theta_1 in range(k_4 + 1):
            for theta_2 in range(k_4 + 1):
                for theta_3 in range(k_2 + 1):
                    f = F(theta_0, theta_1, theta_2, theta_3)
                    if f.sum_of_absolute_values() != 0 and theta_1 + theta_0 + theta_3 + theta_2 == 0:
                            print 4 * "\n"
                            print "3 OJOJOJOJJOOJJOJJ!!!!!!!!!!"
                            print k_0, k_1, k_2, k_3, k_4
                            print theta_0, theta_1, theta_2, theta_3

                    if f.sum_of_absolute_values() == 0 and theta_1 + theta_0 + theta_3 + theta_2 != 0:
                        # print "HURA"
                        # print k_0, k_1, k_2, k_3
                        # print theta_0, theta_1
                        if k_2 != k_3 or theta_0 != theta_1:
                            print 4 * "\n"
                            print "3 SUPER!!!!!!!!!!"
                            print k_0, k_1, k_2, k_3, k_4
                            print theta_0, theta_1, theta_2, theta_3

def tmp(limit=None):
    if limit is None:
        limit = 10
    for k_0 in range(1, limit):
        for k_1 in range(1, limit):
            for k_2 in range(1, limit):
                for k_3 in range(1, limit):
                    first_sum(k_0, k_1, k_2, k_3)
                    for k_4 in range(1, limit):
                        second_sum(k_0, k_1, k_2, k_3, k_4)
