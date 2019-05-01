#Debugging and interacting (put anywhere in code to stop execution for querying structures) (Type: exit() to continue from that point)
#import code
#code.interact(local=locals())

## Testing ##

import knapsack as ns

def size_con(list, size, args):
    if size < 1:
        return False
    else:
        return True

def size_con2(list, size, args):
    if size > 1:
        return False
    else:
        return True

def max_obj(list, size, args):
    total = 0
    for obj in list:
        total += obj.get_value()
    return total

def simple_max(list, size, args):
    total = 0
    for val in list:
        total += val
    return total

def min_obj(list, size, args):
    total = 0
    for obj in list:
        total -= obj.get_value()
    return total

#NOTE: args come in as Tuple, for this case, we only want first tuple item and that is a dictionary
def calories_con(list, size, args):
    cal = 0
    for obj in list:
        if obj.get_name() in args[0]:
            cal += args[0][obj.get_name()]

    if (cal < 750):
        return True
    else:
        return False

def score(list, size, args):
    #args is a tuple with dictionary variables F, ps1 [a], ps2 [b], ps3 [c], ps4 [d], and ps5 [e]
    sum1 = 0
    sum2 = 0
    for obj in list:
        sum2 += obj.get_value()
        if obj.get_name() in args[0]:
            sum1 += args[0][obj.get_name()]*obj.get_value()
    total = 0
    try:
        total = (60.0-sum2)*args[0]["F"] + sum1
    except:
        total = sum1
    return total

def score_con(list, size, args):
    sum = 0
    for obj in list:
        sum += obj.get_value()
    if sum >= 20:
        return True
    else:
        return False


print('\nTest 01')
print('-------')
a = ns.Item()
a.set_name("a")
a.set_value(1)
b = ns.Item()
b.set_name("b")
b.set_value(2)
c = ns.Item()
c.set_name("c")
c.set_value(0)

list = [a,b,c]
for obj in list:
    print(obj)

prob = ns.ZeroOneKnapsack()
prob.register_constraints(size_con)
prob.register_objective_func(min_obj)
prob.exhaustive_search(False)
(val, new_list, status) = prob.Optimize(list)

print('\nAfter Optimization: Test 01 Results')
print(val)
print(status)
print('\nTaken')
for obj in new_list:
    print(obj)


### New Test ###
print('\nTest 02')
print('-------')
foods = []
for i in range(8):
    foods.append(ns.Item())
foods[0].set_name("wine")
foods[0].set_value(89)
foods[1].set_name("beer")
foods[1].set_value(90)
foods[2].set_name("pizza")
foods[2].set_value(95)
foods[3].set_name("burger")
foods[3].set_value(100)
foods[4].set_name("fries")
foods[4].set_value(90)
foods[5].set_name("cola")
foods[5].set_value(79)
foods[6].set_name("apple")
foods[6].set_value(50)
foods[7].set_name("donut")
foods[7].set_value(10)

for obj in foods:
    print(obj)
calories = {"wine":123,"beer":154,"pizza":258,"burger":354,"fries":365,"cola":150,"apple":95,"donut":195}

prob2 = ns.ZeroOneKnapsack()
prob2.register_constraints(calories_con, calories)
prob2.register_objective_func(max_obj)
prob2.exhaustive_search(False)
(val2, new_list2, status2) = prob2.Optimize(foods)

print('\nAfter Optimization: Test 02 Results')
print(val2)
print(status2)
print('\nTaken')
for obj in new_list2:
    print(obj)

### New Test ###
print('\nTest 03')
print('-------')
w = []
for i in range(5):
    w.append(ns.Item())
w[0].set_name("a")
w[0].set_value(10)
w[1].set_name("b")
w[1].set_value(10)
w[2].set_name("c")
w[2].set_value(10)
w[3].set_name("d")
w[3].set_value(10)
w[4].set_name("e")
w[4].set_value(10)
for obj in w:
    print(obj)
grades = {"F":60,"a":10,"b":20,"c":20,"d":10,"e":10}

prob3 = ns.ZeroOneKnapsack()
prob3.register_objective_func(score,grades)
prob3.register_constraints(score_con)
prob3.exhaustive_search(True)
(val3, new_list3, status3) = prob3.Optimize(w)

print('\nAfter Optimization: Test 03 Results')
print(val3)
print(status3)
print('\nTaken')
for obj in new_list3:
    print(obj)

### New Test ###
print('\nTest 04')
print('-------')

vals = [1,2,3]
for obj in vals:
    print(obj)

prob4 = ns.ZeroOneKnapsack()
prob4.register_objective_func(simple_max)
prob4.register_constraints(size_con2)
prob4.exhaustive_search(True)
(val4, new_list4, status4) = prob4.Optimize(vals)

print('\nAfter Optimization: Test 01 Results')
print(val4)
print(status4)
print('\nTaken')
for obj in new_list4:
    print(obj)

# This makes no fucking sense!!!
'''
def test(arg, list = []):
    list.append(arg)
    return list

print(test(1))  # result: [1]        expected: [1]
print(test(2))  # result: [1, 2]     expected: [2]
'''
