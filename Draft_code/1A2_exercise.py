import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')


with open('sydney_rainfall_data.txt') as infile:
    infile.readline()
    infile.readline()
    rainfall = []

    for line in infile:
        lst = line.split()
        lst.pop(0)
        rainfall.append(lst)

    rainfall = [map(float, i) for i in rainfall]


month = 11
amount_partitions = 20
rainfall_by_amount = np.zeros(amount_partitions-1)
monthly_rain = []


for year in rainfall:
    monthly_rain.append(year[month])


partitions = np.linspace(min(map(float, monthly_rain)), max(map(float, monthly_rain)), amount_partitions)


for rain in monthly_rain:
    for i in xrange(0, amount_partitions-1):
        if rain >= partitions[i] and rain <= partitions[i+1]:
            rainfall_by_amount[i] += 1
            break


normalized_rainfall_amount = (rainfall_by_amount - min(rainfall_by_amount)) / (np.sum(rainfall_by_amount))
"""
plt.xticks(xrange(len(normalized_rainfall_amount)))
plt.bar(xrange(len(normalized_rainfall_amount)), normalized_rainfall_amount)
plt.show()
"""


monthly_average = np.zeros(11)
for i in xrange(11):
    summed_rain = 0
    for year in rainfall:
        summed_rain += year[i]
    monthly_average[i] = summed_rain / len(rainfall)

print partitions
