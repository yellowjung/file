import re
import sys

filename = sys.argv[1]
data = open(filename, 'r', errors='replace')
result = [0,0,0,0,0]
count = 0
r_text = 0
error = []
empty = []
while True:
    readTmp = data.readline()
    if not readTmp : break
    Tmp = readTmp.split(',')
    try:
        if '"empty"' in Tmp[2]:
            empty.append(count)
            continue
        entropy = float(Tmp[2].replace("\"", ""))
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

    except:
        error.append(count)
        print(count, Tmp)

    count += 1

print("Error text count: ", len(error))
print("Empty text count: ", len(empty))

for res in result:
    r_text += res
print("Correct text count: ", r_text)
print("Error ratio : ",len(error)/(count+len(error)+len(empty))*100,"%")
print("Empty ratio : ",len(empty)/(count+len(error)+len(empty))*100,"%")
print("Result:")
for res in result:
    print(res)
print("Total count : ", count + len(error) + len(empty))
data.close()
