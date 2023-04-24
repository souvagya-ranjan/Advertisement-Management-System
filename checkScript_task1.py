from collections import defaultdict
def comparefiles(file1, file2, k1, k2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        for i in range(k1, k2):
            val1 = int(f1.readline().strip())
            val2 = int(f2.readline().strip())
            if val1!=val2:
                print("output for k does not match", i, val1, val2)
                return False
            if val1==0:
                continue  
            size1 = int(f1.readline().strip())
            size2 = int(f2.readline().strip())
            if size1!=size2:
                print("output for size does not match for k ", i, size1, size2)
                return False
            arr1 = []
            arr2 = []
            for j in range(size1):
                ele1 = f1.readline().strip()
                ele2 = f2.readline().strip()
                arr1.append(list(map(int, ele1.split())))
                arr2.append(list(map(int, ele2.split())))
            arr1.sort(key = lambda x: x[0])
            arr2.sort(key = lambda x: x[0])
            for j in range(size1):
                if arr1[j]!=arr2[j]:
                    print("output for array does not match for k", i, j,arr1[j], arr2[j])
                    return False
    return True

if __name__ == "__main__":
    file1 = "A3_test/test3/task2_output3_verbose.txt"
    file2 = "Output/task2_output3.txt"
    k1 = 1
    k2 = 5
    print(comparefiles(file1, file2, k1, k2))
        