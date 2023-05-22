import math as m

def Stability_Check_1():

    s = 2.0
    d = m.sqrt(2)
    for i in range(1, 200):
        print()
        d *= m.sqrt(2)
        print(f"i = {i}")
        print(f"d = {d}")
        print(f"s = {s}")
        print(f"(s-d)   = {(s-d)}")
        print(f"(s-d)/s = {(s-d)/s}")
        d *= m.sqrt(2)
        s *= 2.0

    return

def Stability_Check_2():

    s = 3.0
    d = m.sqrt(3)
    for i in range(1, 200):
        print()
        d *= m.sqrt(3)
        print(f"i = {i}")
        print(f"d = {d}")
        print(f"s = {s}")
        print(f"(s-d)   = {(s-d)}")
        print(f"(s-d)/s = {(s-d)/s}")
        d *= m.sqrt(3)
        s *= 3.0

    return

import gmpy2 as mp

def Stability_Check_MultiPrecision(floatnumber, precision, numberOfIterations):

    mp.set_context(mp.context())
    mp.get_context().precision = precision

    s = mp.mpfr(floatnumber)            # ground truth
    d = mp.sqrt(mp.mpfr(floatnumber))

    for i in range(1, numberOfIterations):
        d = mp.mul(d, mp.sqrt(mp.mpfr(floatnumber)))
        d = mp.mul(d, mp.sqrt(mp.mpfr(floatnumber)))
        s = mp.mul(s, mp.mpfr(floatnumber))

    d = mp.mul(d, mp.sqrt(mp.mpfr(floatnumber)))
    ratio = mp.div(mp.sub(s, d), s)

    return ratio

def Stability_Check_MultiPrecision_Print(floatnumber, precision, printon):

    mp.set_context(mp.context())
    mp.get_context().precision = precision

    s = mp.mpfr(floatnumber)

    if printon:
        print(s.precision)

    d = mp.sqrt(mp.mpfr(floatnumber))

    for i in range(1, 100):
        if printon:
            print()
        d = mp.mul(d, mp.sqrt(mp.mpfr(floatnumber)))
        if printon:
            print("i = ", i)
            print("d = ", d)
            print("s = ", s)

        diff_float = mp.sub(s, d)
        if printon:
            print("(s-d)   = ", diff_float)

        ratio = mp.div(diff_float, s)
        if printon:
            print("(s-d)/s = ", ratio)

        d = mp.mul(d, mp.sqrt(mp.mpfr(floatnumber)))
        s = mp.mul(s, mp.mpfr(floatnumber))

    if True: #printon:
        d = mp.mul(d, mp.sqrt(mp.mpfr(floatnumber)))
        print("s = ", s)
        print("d = ", d)

    return ratio


import pandas as pd

def Test_Stability():

    floatnumbers = [3.7]
    iterations   = [300]
    precisions = []

    for precision in range(1, 100, 1):
        precisions.append(precision)

    floatnumber_list = []
    iteration_list = []
    precision_list = []
    ratio_list = []

    for floatnumber in floatnumbers:
        for iteration in iterations:
            for precision in precisions:

                #ratio = mp.log10(Stability_Check_MultiPrecision(floatnumber, precision, iteration))
                ratio = Stability_Check_MultiPrecision(floatnumber, precision, iteration)

                floatnumber_list.append(floatnumber)
                precision_list.append(precision)
                iteration_list.append(iteration)
                ratio_list.append(ratio)

    rows = ({'floatnumber': floatnumber_list, 'precision': precision_list, 'iteration': iteration_list, 'ratio': ratio_list})
    df = pd.DataFrame(rows)

    return df

#----------------------------------------

df = Test_Stability()

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

ax.plot(df['precision'], df['ratio'])
ax.set_title('Ratio vs Precision with iter = 300')
ax.set_xlabel('precision [number of digits]')
ax.set_ylabel('ratio log10[(s-d)/s]')

fig.savefig('ratio_precision - iter 300 pre 1-100.pdf')
fig.show()

df.to_csv('ratio_precision - iter 300 pre 1-100.csv')

#----------------------------------------

iterations = 300

start =    1
end   =  100
step  =    1

print("Floatnumber = ", 1.0)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(1.0, i, iterations))
print()

print("Floatnumber = ", 2.0)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(2.0, i, iterations))
print()

print("Floatnumber = ", 4.0)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(4.0, i, iterations))
print()

print("Floatnumber = ", 1.1)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(1.1, i, iterations))
print()

print("Floatnumber = ", 2.2)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(2.2, i, iterations))
print()

print("Floatnumber = ", 4.4)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(4.4, i, iterations))
print()

print("Floatnumber = ", 1.7)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(1.7, i, iterations))
print()

print("Floatnumber = ", 5.1)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(5.1, i, iterations))
print()

print("Floatnumber = ", 9.1)
for i in range(start, end, step):
    print("Precision, (d-s)/s = ", i, Stability_Check_MultiPrecision(9.1, i, iterations))
print()

