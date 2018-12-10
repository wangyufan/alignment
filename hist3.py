from numpy import array
from numpy.random import normal
from matplotlib import pyplot
import matplotlib.pyplot as plt

def get_data(filename):
    n=0
    data=[]
    readfile=open(filename)
    lines=readfile.readlines()
    for line in lines:
        if n>0 and n<4000 :
            line=line.split()[6].split(")")[0]
            data.append(str(line))
        n+=1
    res = (data.count("'Active'"), data.count("'Moderate'"), data.count("'Inactive'"))
    print res
    return res

data_595 = get_data("all_new_anticc595.txt")
data_763 = get_data("all_new_anticc763.txt")
data_1452 = get_data("all_new_anticc1452.txt")

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/4., 1.01*height, '%s' % float(height))

if __name__ == '__main__': 

    pyplot.xlabel('Activity Type')
    pyplot.ylabel('Frequency')
    pyplot.title('Activity Frequency of Top40 Molecules selected')
    name = ["'Active'" "'Moderate'" "'Inactive'"]
    total_width, n = 0.9, 3 
    width = total_width / n
    x=[0,1,2]

    bar_595 = plt.bar(x, data_595, width=width, fc = 'y', label= 'alldata_595') 
    for i in range(len(x)):  
        x[i] = x[i] + width 
    bar_763 = plt.bar(x, data_763, width=width, fc = 'g', label= 'alldata_763') 
    for i in range(len(x)):  
        x[i] = x[i] + width 
    bar_1452 = plt.bar(x, data_1452, width=width, fc = 'b', label= 'alldata_1452') 

    plt.xticks((0.3,1.3,2.3),("Active", "Moderate", "Inactive"))

    autolabel(bar_595)
    autolabel(bar_763)
    autolabel(bar_1452)
    plt.legend()
    pyplot.show()