import psutil

for i in range(0,10000):
    a = 100000000
    a = 0
    if i % 100 == 0:
        print(psutil.virtual_memory())
