import sympy as sp
from sympy.abc import x
import tkinter as tk
from tkinter import messagebox

def string_to_list(string):
    """

    Parameters
    ----------
    string : str
        A string

    Returns
    -------
    array : list (str)
        Characters stored in a array

    """
    array = string.split()
    for i in range(0, len(array)):
        array[i] = float(array[i])
    return(array)

def absolute(num):
    """
    
    Parameters
    ----------
    num : float
        A floating point number

    Returns
    -------
    The absolute value of num

    """
    if num < 0:
        return(-num)
    else:
        return(num)
    
def linspace(a, b, n):
    """

    Parameters
    ----------
    a : float
        The beginning of the interval
    b : float
        The end of the interval
    n : int
        The number of subintervals desired

    Returns
    -------
    linspace_array : list (float)
        An array of numbers from a to b such that there are n subintervals

    """
    distance = (b - a) / n
    linspace_array = [a]
    for i in range(1, n + 1):
        linspace_array += [a + i * distance]
    return(linspace_array)

def bisection_method(f, a, b, tolerance, max_num):
    """

    Parameters
    ----------
    f : sp.expr
        A sympy expression
    a : float
        The beginning of the interval of interest
    b : float
        The end of the interval of interest
    tolerance : float
        Maximum tolerated error
    max_num : int
        The maximum number of iterations desired

    Returns
    -------
    A string message containing either a warning or an approximate root of
    the function

    """
    i = 0
    fa = f.evalf(subs = {x:a})
    fb = f.evalf(subs = {x:b})
    if fa * fb > 0:
        return("Both f(a) and f(b) are the same sign. Please try again.")
    while i < max_num:
        p = (a + b) / 2
        fp = f.evalf(subs = {x:p})
        if absolute(fp) < tolerance:
            return(str(p) + " is an approximate root of the function.")
        fa = f.evalf(subs = {x:a})
        if fa * fp > 0:
            a = p
        else:
            b = p
        i += 1
    return("The Bisection Method could not produce a root in under " +
           str(max_num) + " iterations for the tolerance level given.")

def fixed_point_iteration(f, p, tolerance, max_num):
    """

    Parameters
    ----------
    f : sp.expr
        A sympy expression
    p : float
        An initial guess for the fixed point
    tolerance : float
        Maximum tolerated error
    max_num : int
        The maximum number of iterations desired

    Returns
    -------
    A string message containing either a warning or an approximate fixed
    point of the function

    """
    i = 0
    while i < max_num:
        fp = f.evalf(subs = {x:p})
        if absolute(fp - p) < tolerance:
            return(str(p) + " is an approximate fixed point of the function.")
        p = fp
        i += 1
    return("Fixed Point Iteration could not produce a fixed point in under " +
           str(max_num) + " iterations for the tolerance level given.")

def lagrange_polynomial(f_array, x_array):
    """

    Parameters
    ----------
    f_array : list (float)
        An array of f(x) values
    x_array : list (float)
        An array of x values

    Returns
    -------
    A string message containing either a warning or an interpolated polynomial

    """
    if len(f_array) != len(x_array):
        return("The length of both arrays must be equivalent to one another.")
    poly_num_array = [""] * len(f_array)
    poly_denom_array = [1] * len(f_array)
    poly = ""
    for i in range(0, len(f_array)):
        for j in range(0, len(f_array)):
            if j != i:
                poly_num_array[i] += "(x - " + str(x_array[j]) + ")"
                poly_denom_array[i] *= (x_array[i] - x_array[j])
    for k in range(0, len(f_array)):
        poly += poly_num_array[k]
        poly += str(f_array[k]) + "/" + str(poly_denom_array[k])
        if k != len(f_array) - 1:
            poly += " + "
    return(poly)

def three_point_midpoint_derivative(f_array, x_array):
    """

    Parameters
    ----------
    f_array : list (float)
        An array of f(x) values
    x_array : list (float)
        An array of x values

    Returns
    -------
    A string message containing either a warning or an approximation of the
    first derivative at the given point

    """
    if (len(f_array) != 3) or (len(x_array) != 3):
        return("The length of both arrays must be 3.")
    elif (x_array[2] - x_array[1]) != (x_array[1] - x_array[0]):
        return("The spacing of the x values must be consistent.")
    else:
        h = x_array[2] - x_array[1]
        f_prime = (1 / (h * 2)) * (f_array[2] - f_array[0])
        return("The derivative of the function at f is approximately " +
               str(f_prime) + " at " + str(x_array[1]) + ".")
    
def three_point_endpoint_derivative(f_array, x_array):
    """

    Parameters
    ----------
    f_array : list (float)
        An array of f(x) values
    x_array : list (float)
        An array of x values

    Returns
    -------
    A string message containing either a warning or an approximation of the
    first derivative at the given point

    """
    if (len(f_array) != 3) or (len(x_array) != 3):
        return("The length of both arrays must be 3.")
    elif (x_array[2] - x_array[1]) != (x_array[1] - x_array[0]):
        return("The spacing of the x values must be consistent.")
    else:
        h = x_array[2] - x_array[1]
        f_prime = (1 / (h * 2)) * (-3 * f_array[0] + 4 * f_array[1]
                                   - f_array[2])
        return("The derivative of the function at f is approximately " +
               str(f_prime) + " at " + str(x_array[0]) + ".")
    
def composite_simpsons_rule(f, a, b, n):
    """

    Parameters
    ----------
    f : sp.expr
        A sympy expression
    a : float
        The beginning of the interval of integration
    b : float
        The end of the interval of integration
    n : int
        The number of subintervals desired

    Returns
    -------
    A string message containing either a warning or an approximation of the
    integral over the given interval

    """
    if n % 2 != 0:
        return("There must be an even number of subintervals.")
    else:
        x_array = linspace(a, b, n)
        first_sum = 0
        for i in range(1, int(n / 2)):
            first_sum += 2 * f.evalf(subs = {x:x_array[2 * i]})
        second_sum = 0
        for i in range(1, int(n / 2 + 1)):
            second_sum += 4 * f.evalf(subs = {x:x_array[2 * i - 1]})
        h = (a - b) / n
        f0 = f.evalf(subs = {x:a})
        fn = f.evalf(subs = {x:b})
        integral = (h / 3) * (f0 + first_sum + second_sum + fn)
        return("The integral of f over the region is " +
               str(-integral) + ".")
    
def composite_trapezoidal_rule(f, a, b, n):
    """

    Parameters
    ----------
    f : sp.expr
        A sympy expression
    a : float
        The beginning of the interval of integration
    b : float
        The end of the interval of integration
    n : int
        The number of subintervals desired

    Returns
    -------
    A string message containing either a warning or an approximation of the
    integral over the given interval

    """
    x_array = linspace(a, b, n)
    first_sum = 0
    for i in range(1, int(n)):
        first_sum += 2 * f.evalf(subs = {x:x_array[i]})
    h = (a - b) / n
    f0 = f.evalf(subs = {x:a})
    fn = f.evalf(subs = {x:b})
    integral = (h / 2) * (f0 + first_sum + fn)
    return("The integral of f over the region is " +
               str(-integral) + ".")

# Designed for up to 5 arguments, can be modified to add more
# If adding another numerical method, add instructions within the
# return_answer() function as to how to treat the string user input

def second_window(window, method, inputs, method_name):
    window_2 = tk.Toplevel(window)
    window_2.title(method_name)
    window_2.geometry('500x500')
    window_2.grab_set()
    args = [""] * 5
    strings = []
    for i in range(0, 5):
        strings += [tk.StringVar()]
    def retrieve_arg1():
        str_input = strings[0].get()
        args[0] = str_input
    def retrieve_arg2():
        str_input = strings[1].get()
        args[1] = str_input
    def retrieve_arg3():
        str_input = strings[2].get()
        args[2] = str_input
    def retrieve_arg4():
        str_input = strings[3].get()
        args[3] = str_input
    def retrieve_arg5():
        str_input = strings[4].get()
        args[4] = str_input
    fun_list = [retrieve_arg1, retrieve_arg2, retrieve_arg3,
                retrieve_arg4, retrieve_arg5]
    def return_answer():
        if method_name == "Bisection Method":
            answer = bisection_method(sp.sympify(args[0]),
                                      float(args[1]), float(args[2]),
                                      float(args[3]), int(args[4]))
            messagebox.showinfo("Answer", answer)
        elif method_name == "Fixed Point Iteration":
            answer = fixed_point_iteration(sp.sympify(args[0]),
                                           float(args[1]),
                                           float(args[2]), int(args[3]))
            messagebox.showinfo("Answer", answer)
        elif method_name == "Lagrange Polynomial":
            f_array = string_to_list(args[0])
            x_array = string_to_list(args[1])
            answer = lagrange_polynomial(f_array, x_array)
            messagebox.showinfo("Answer", answer)
        elif method_name == "Three Point Midpoint (Differentiation)":
            f_array = string_to_list(args[0])
            x_array = string_to_list(args[1])
            answer = three_point_midpoint_derivative(f_array, x_array)
            messagebox.showinfo("Answer", answer)
        elif method_name == "Three Point Endpoint (Differentiation)":
            f_array = string_to_list(args[0])
            x_array = string_to_list(args[1])
            answer = three_point_endpoint_derivative(f_array, x_array)
            messagebox.showinfo("Answer", answer)
        elif method_name == "Composite Simpson's Rule":
            answer = composite_simpsons_rule(sp.sympify(args[0]),
                                             float(args[1]), float(args[2]),
                                             int(args[3]))
            messagebox.showinfo("Answer", answer)
        elif method == "Composite Trapezoidal Rule":
            answer = composite_trapezoidal_rule(sp.sympify(args[0]),
                                                float(args[1]),
                                                float(args[2]),
                                                int(args[3]))
            messagebox.showinfo("Answer", answer)
    labels = []
    entries = []
    buttons = []
    for i in range(0, 5):
        labels += [tk.Label(window_2, text = inputs[int(method)][i])]
        labels[i].pack()
        entries += [tk.Entry(window_2, textvariable = strings[i])]
        entries[i].pack()
        buttons += [tk.Button(window_2, text = "Enter", command = fun_list[i])]
        buttons[i].pack()
    button = tk.Button(window_2, text = "Compute", command = return_answer)
    button.pack()
    window_2.mainloop()

def initial_window(methods):
    window = tk.Tk()
    window.title("Numerical Method Calculator")
    window.geometry("500x250")
    str_var = tk.StringVar()
    Radiobuttons = []
    def retrieve_input():
        str_input = str_var.get()
        second_window(window, str_input, inputs, methods[int(str_input)])
    label = tk.Label(window, text = "Please select a Numerical Method.")
    label.pack()
    for i in range(0, len(methods)):
        Radiobuttons += [tk.Radiobutton(window, text = methods[i],
                                        variable = str_var, value = i)]
        Radiobuttons[i].pack()
    button = tk.Button(window, text = "Continue", command = retrieve_input)
    button.pack()
    window.mainloop()
    
methods = ["Bisection Method", "Fixed Point Iteration",
           "Lagrange Polynomial", "Three Point Midpoint (Differentiation)",
           "Three Point Endpoint (Differentiation)",
           "Composite Simpson's Rule", "Composite Trapezoidal Rule"]
    
inputs = [["Enter the function f(x), e.g. sin(x) + 3*x**2 + 13 / (x - 2).",
           "Set a left bound for the search.",
           "Set a right bound for the search.",
           "Set an error bound.",
           "Set a maximum number of iterations."]]
inputs += [["Enter the function f(x), e.g. sin(x) + 3*x**2 + 13 / (x - 2).",
            "Set an initial guess for a fixed point.",
            "Set an error bound.",
            "Set a maximum number of iterations.",
            "Leave empty."]]
inputs += [["Enter an array of f(x) values, e.g. 3 3.2 4.",
            "Enter an array of x values, e.g. 1 2 3.6.",
            "Leave empty.",
            "Leave empty.",
            "Leave empty."]]
inputs += [["Enter an array of f(x) values, e.g. 1 3.2 4.",
            "Enter an array of evenly spaced x values, e.g. 1 2 3 or 3 2 1.",
            "Leave empty.",
            "Leave empty.",
            "Leave empty."]] * 2
inputs += [["Enter the function f(x), e.g. sin(x) + 3*x**2 + 13 / (x - 2).",
            "Set an a for the interval of integration.",
            "Set a b for the interval of integration.",
            "Set a number of intervals (even for Simpson's).",
            "Leave empty."]] * 2

initial_window(methods)