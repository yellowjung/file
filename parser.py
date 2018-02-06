import re
import sys

filename = sys.argv[1]
data = open(filename)
result = [0,0,0,0,0]
count = 0
r_text = 0
error = []
while True:
    readTmp = data.readline()
    Tmp = readTmp.split(',')
    try:
        entropy = float(Tmp[2])
    except:
        error.append(count)
    if entropy < 5:
        result[0] += 1
    elif entropy >= 5 and entropy < 6:
        result[1] += 1
    elif entropy >= 6 and entropy < 7:
        result[2] += 1
    elif entropy >= 7 and entropy < 8:
        result[3] += 1
    elif entropy >= 8:
        result[4] += 1
    if not readTmp : break
    count += 1

print("Error text count: ", len(error))

for res in result:
    r_text += res
print("Correct text count: ", r_text)
print("ratio : ",len(error)/count*100,"%")
print("Result:")
for res in result:
    print(res)
data.close()
